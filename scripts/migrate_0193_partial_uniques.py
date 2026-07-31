#!/usr/bin/env python3
r"""Apply 0193 link-table partial unique indexes table-by-table via psycopg2.

For each table, independently:
  1. Drop legacy unique_together constraints (short AccessExclusiveLock)
  2. Deduplicate NULL-feature rows (DELETE … USING)
  3. CREATE UNIQUE INDEX CONCURRENTLY for both partial uniques

Uses autocommit so CONCURRENTLY is legal and locks are not held across tables.

Examples:
  python scripts/migrate_0193_partial_uniques.py \\
      --dsn "postgresql://user:pass@host:5432/dbname"

  python scripts/migrate_0193_partial_uniques.py --dsn "$DATABASE_URL" --table lamindb_artifactulabel

  # Status only (no changes):
  python scripts/migrate_0193_partial_uniques.py --dsn "$DATABASE_URL" --stat

  # After SQL succeeds, record Django migration as applied:
  python scripts/migrate_0193_partial_uniques.py --dsn "$DATABASE_URL" --fake-django
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass

import psycopg2
from psycopg2 import sql

MIGRATION_APP = "lamindb"
MIGRATION_NAME = "0193_v2_9_part_2"


@dataclass(frozen=True)
class LinkTable:
    table: str
    pair: tuple[str, str]
    triple: tuple[str, ...]
    name: str
    name_null: str


LINK_TABLES: list[LinkTable] = [
    LinkTable(
        "lamindb_artifactartifact",
        ("artifact_id", "value_id"),
        ("artifact_id", "value_id", "feature_id"),
        "unique_artifactartifact",
        "unique_artifactartifact_null_feature",
    ),
    LinkTable(
        "lamindb_artifactproject",
        ("artifact_id", "project_id"),
        ("artifact_id", "project_id", "feature_id"),
        "unique_artifactproject",
        "unique_artifactproject_null_feature",
    ),
    LinkTable(
        "lamindb_artifactrecord",
        ("artifact_id", "record_id"),
        ("artifact_id", "record_id", "feature_id"),
        "unique_artifactrecord",
        "unique_artifactrecord_null_feature",
    ),
    LinkTable(
        "lamindb_artifactreference",
        ("artifact_id", "reference_id"),
        ("artifact_id", "reference_id", "feature_id"),
        "unique_artifactreference",
        "unique_artifactreference_null_feature",
    ),
    LinkTable(
        "lamindb_artifactrun",
        ("artifact_id", "run_id"),
        ("artifact_id", "run_id", "feature_id"),
        "unique_artifactrun",
        "unique_artifactrun_null_feature",
    ),
    LinkTable(
        "lamindb_artifactulabel",
        ("artifact_id", "ulabel_id"),
        ("artifact_id", "ulabel_id", "feature_id"),
        "unique_artifactulabel",
        "unique_artifactulabel_null_feature",
    ),
    LinkTable(
        "lamindb_artifactuser",
        ("artifact_id", "user_id"),
        ("artifact_id", "user_id", "feature_id"),
        "unique_artifactuser",
        "unique_artifactuser_null_feature",
    ),
    LinkTable(
        "lamindb_collectionrecord",
        ("collection_id", "record_id"),
        ("collection_id", "record_id", "feature_id"),
        "unique_collectionrecord",
        "unique_collectionrecord_null_feature",
    ),
    LinkTable(
        "lamindb_projectrecord",
        ("project_id", "record_id"),
        ("project_id", "feature_id", "record_id"),
        "unique_projectrecord",
        "unique_projectrecord_null_feature",
    ),
    LinkTable(
        "lamindb_referencerecord",
        ("reference_id", "record_id"),
        ("reference_id", "feature_id", "record_id"),
        "unique_referencerecord",
        "unique_referencerecord_null_feature",
    ),
    LinkTable(
        "lamindb_runartifact",
        ("run_id", "artifact_id"),
        ("run_id", "artifact_id", "feature_id"),
        "unique_runartifact",
        "unique_runartifact_null_feature",
    ),
    LinkTable(
        "lamindb_runrecord",
        ("run_id", "record_id"),
        ("run_id", "record_id", "feature_id"),
        "unique_runrecord",
        "unique_runrecord_null_feature",
    ),
    LinkTable(
        "lamindb_transformrecord",
        ("transform_id", "record_id"),
        ("transform_id", "record_id", "feature_id"),
        "unique_transformrecord",
        "unique_transformrecord_null_feature",
    ),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def quote_ident(name: str) -> sql.Identifier:
    return sql.Identifier(name)


def drop_legacy_uniques(cur, table: str, *, dry_run: bool) -> list[str]:
    cur.execute(
        """
        SELECT c.conname
        FROM pg_constraint c
        JOIN pg_class t ON c.conrelid = t.oid
        JOIN pg_namespace n ON t.relnamespace = n.oid
        WHERE n.nspname = current_schema()
          AND t.relname = %s
          AND c.contype = 'u'
        """,
        (table,),
    )
    names = [row[0] for row in cur.fetchall()]
    for conname in names:
        stmt = sql.SQL("ALTER TABLE {} DROP CONSTRAINT IF EXISTS {}").format(
            quote_ident(table), quote_ident(conname)
        )
        log(f"  DROP CONSTRAINT {conname}")
        if not dry_run:
            cur.execute(stmt)
    return names


def dedupe_null_feature(cur, spec: LinkTable, *, dry_run: bool) -> int:
    a, b = spec.pair
    stmt = sql.SQL(
        """
        DELETE FROM {table} AS keep
        USING {table} AS dup
        WHERE keep.feature_id IS NULL
          AND dup.feature_id IS NULL
          AND keep.{a} = dup.{a}
          AND keep.{b} = dup.{b}
          AND keep.id > dup.id
        """
    ).format(
        table=quote_ident(spec.table),
        a=quote_ident(a),
        b=quote_ident(b),
    )
    if dry_run:
        # Count only for dry-run visibility.
        count_stmt = sql.SQL(
            """
            SELECT COUNT(*) FROM {table} AS keep
            JOIN {table} AS dup
              ON keep.feature_id IS NULL
             AND dup.feature_id IS NULL
             AND keep.{a} = dup.{a}
             AND keep.{b} = dup.{b}
             AND keep.id > dup.id
            """
        ).format(
            table=quote_ident(spec.table),
            a=quote_ident(a),
            b=quote_ident(b),
        )
        cur.execute(count_stmt)
        return int(cur.fetchone()[0])

    cur.execute(stmt)
    return cur.rowcount


def index_exists(cur, name: str) -> bool:
    cur.execute("SELECT 1 FROM pg_class WHERE relname = %s AND relkind = 'i'", (name,))
    return cur.fetchone() is not None


def index_is_valid(cur, name: str) -> bool:
    status = index_status(cur, name)
    return status == "valid"


def index_status(cur, name: str) -> str:
    """Return 'valid', 'invalid', or 'missing'."""
    cur.execute(
        """
        SELECT i.indisvalid
        FROM pg_class c
        JOIN pg_index i ON i.indexrelid = c.oid
        WHERE c.relname = %s AND c.relkind = 'i'
        """,
        (name,),
    )
    row = cur.fetchone()
    if row is None:
        return "missing"
    return "valid" if row[0] else "invalid"


def legacy_unique_constraints(cur, table: str) -> list[str]:
    cur.execute(
        """
        SELECT c.conname
        FROM pg_constraint c
        JOIN pg_class t ON c.conrelid = t.oid
        JOIN pg_namespace n ON t.relnamespace = n.oid
        WHERE n.nspname = current_schema()
          AND t.relname = %s
          AND c.contype = 'u'
        """,
        (table,),
    )
    return [row[0] for row in cur.fetchall()]


def count_null_feature_dupes(cur, spec: LinkTable) -> int:
    a, b = spec.pair
    stmt = sql.SQL(
        """
        SELECT COUNT(*) FROM {table} AS keep
        JOIN {table} AS dup
          ON keep.feature_id IS NULL
         AND dup.feature_id IS NULL
         AND keep.{a} = dup.{a}
         AND keep.{b} = dup.{b}
         AND keep.id > dup.id
        """
    ).format(
        table=quote_ident(spec.table),
        a=quote_ident(a),
        b=quote_ident(b),
    )
    cur.execute(stmt)
    return int(cur.fetchone()[0])


def django_migration_recorded(cur) -> bool:
    cur.execute(
        """
        SELECT 1 FROM django_migrations
        WHERE app = %s AND name = %s
        """,
        (MIGRATION_APP, MIGRATION_NAME),
    )
    return cur.fetchone() is not None


def print_stats(conn, tables: list[LinkTable]) -> None:
    """Read-only progress report for 0193 partial uniques."""
    n_indexes = 0
    n_valid = 0
    n_invalid = 0
    n_missing = 0
    n_legacy = 0
    n_dupes = 0
    n_tables_done = 0

    log("\n=== status ===")
    with conn.cursor() as cur:
        recorded = django_migration_recorded(cur)
        log(
            f"django_migrations {MIGRATION_APP}.{MIGRATION_NAME}: "
            f"{'recorded' if recorded else 'not recorded'}"
        )
        log(
            f"{'table':<32} {'legacy':>6} {'dupes':>7} "
            f"{'idx_feature':<10} {'idx_null':<10} {'ready':>5}"
        )
        log("-" * 80)

        for spec in tables:
            legacy = legacy_unique_constraints(cur, spec.table)
            dupes = count_null_feature_dupes(cur, spec)
            st_feat = index_status(cur, spec.name)
            st_null = index_status(cur, spec.name_null)

            for st in (st_feat, st_null):
                n_indexes += 1
                if st == "valid":
                    n_valid += 1
                elif st == "invalid":
                    n_invalid += 1
                else:
                    n_missing += 1

            n_legacy += len(legacy)
            n_dupes += dupes
            ready = (
                not legacy and dupes == 0 and st_feat == "valid" and st_null == "valid"
            )
            if ready:
                n_tables_done += 1

            log(
                f"{spec.table:<32} {len(legacy):>6} {dupes:>7} "
                f"{st_feat:<10} {st_null:<10} {'yes' if ready else 'no':>5}"
            )

    log("-" * 80)
    log(
        f"tables ready: {n_tables_done}/{len(tables)}  |  "
        f"indexes valid/invalid/missing: {n_valid}/{n_invalid}/{n_missing} "
        f"(of {n_indexes})  |  "
        f"legacy constraints: {n_legacy}  |  "
        f"pending NULL-feature dupes: {n_dupes}"
    )


def create_unique_index_concurrently(
    conn,
    *,
    name: str,
    table: str,
    columns: tuple[str, ...],
    where_sql: str,
    dry_run: bool,
) -> None:
    old = conn.autocommit
    conn.autocommit = True
    try:
        with conn.cursor() as cur:
            if index_exists(cur, name) and index_is_valid(cur, name):
                log(f"  index {name}: already valid, skip create")
                return
            if index_exists(cur, name):
                log(f"  index {name}: INVALID leftover, dropping first")
                if not dry_run:
                    cur.execute(
                        sql.SQL("DROP INDEX CONCURRENTLY IF EXISTS {}").format(
                            quote_ident(name)
                        )
                    )

            cols = sql.SQL(", ").join(quote_ident(c) for c in columns)
            stmt = sql.SQL(
                "CREATE UNIQUE INDEX CONCURRENTLY {} ON {} ({}) WHERE {}"
            ).format(
                quote_ident(name),
                quote_ident(table),
                cols,
                sql.SQL(where_sql),
            )
            log(f"  CREATE UNIQUE INDEX CONCURRENTLY {name}")
            if not dry_run:
                cur.execute(stmt)
    finally:
        conn.autocommit = old


def migrate_table(conn, spec: LinkTable, *, dry_run: bool) -> None:
    started = time.perf_counter()
    log(f"\n=== {spec.table} ===")

    # Steps that need a normal transaction: drop constraint + dedupe.
    conn.autocommit = False
    try:
        with conn.cursor() as cur:
            dropped = drop_legacy_uniques(cur, spec.table, dry_run=dry_run)
            if not dropped:
                log("  no legacy unique constraints")

            deleted = dedupe_null_feature(cur, spec, dry_run=dry_run)
            log(f"  dedupe NULL-feature rows: {deleted}")

        if dry_run:
            conn.rollback()
        else:
            conn.commit()
    except Exception:
        conn.rollback()
        raise

    # Index DDL must be outside a transaction.
    create_unique_index_concurrently(
        conn,
        name=spec.name,
        table=spec.table,
        columns=spec.triple,
        where_sql='"feature_id" IS NOT NULL',
        dry_run=dry_run,
    )
    create_unique_index_concurrently(
        conn,
        name=spec.name_null,
        table=spec.table,
        columns=spec.pair,
        where_sql='"feature_id" IS NULL',
        dry_run=dry_run,
    )

    elapsed = time.perf_counter() - started
    log(f"  done in {elapsed:.1f}s")


def fake_django_migration(conn, *, dry_run: bool) -> None:
    log(f"\n=== fake Django migration {MIGRATION_APP}.{MIGRATION_NAME} ===")
    conn.autocommit = False
    try:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT 1 FROM django_migrations
                WHERE app = %s AND name = %s
                """,
                (MIGRATION_APP, MIGRATION_NAME),
            )
            if cur.fetchone():
                log("  already recorded, skip")
                conn.rollback()
                return
            log("  INSERT INTO django_migrations ...")
            if not dry_run:
                cur.execute(
                    """
                    INSERT INTO django_migrations (app, name, applied)
                    VALUES (%s, %s, NOW())
                    """,
                    (MIGRATION_APP, MIGRATION_NAME),
                )
                conn.commit()
            else:
                conn.rollback()
    except Exception:
        conn.rollback()
        raise


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--dsn",
        default=os.environ.get("DATABASE_URL")
        or os.environ.get("LAMINDB_DJANGO_DATABASE_URL"),
        help="Postgres DSN (or set DATABASE_URL / LAMINDB_DJANGO_DATABASE_URL)",
    )
    p.add_argument(
        "--table",
        action="append",
        dest="tables",
        help="Only migrate this table (repeatable). Default: all link tables.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions / counts without changing the database.",
    )
    p.add_argument(
        "--stat",
        action="store_true",
        help="Print per-table status (indexes/legacy/dupes) and exit without changes.",
    )
    p.add_argument(
        "--fake-django",
        action="store_true",
        help=f"Record {MIGRATION_APP}.{MIGRATION_NAME} in django_migrations after SQL.",
    )
    p.add_argument(
        "--lock-timeout",
        default="5s",
        help="Postgres lock_timeout for DROP CONSTRAINT / DELETE (default: 5s).",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.dsn:
        log("error: pass --dsn or set DATABASE_URL")
        return 2

    selected = LINK_TABLES
    if args.tables:
        wanted = set(args.tables)
        selected = [t for t in LINK_TABLES if t.table in wanted]
        missing = wanted - {t.table for t in selected}
        if missing:
            log(f"error: unknown table(s): {sorted(missing)}")
            return 2

    log("connecting…")
    conn = psycopg2.connect(args.dsn)
    try:
        conn.autocommit = False
        with conn.cursor() as cur:
            cur.execute("SELECT current_database(), current_user")
            db, user = cur.fetchone()
            log(f"database={db} user={user}")
            if args.lock_timeout and not args.stat:
                cur.execute(
                    "SELECT set_config('lock_timeout', %s, false)",
                    (args.lock_timeout,),
                )
            conn.commit()

        if args.stat:
            print_stats(conn, selected)
            return 0

        for spec in selected:
            # Retry DROP/DELETE a few times if lock_timeout fires under load.
            attempts = 0
            while True:
                attempts += 1
                try:
                    migrate_table(conn, spec, dry_run=args.dry_run)
                    break
                except psycopg2.errors.LockNotAvailable:
                    conn.rollback()
                    if attempts >= 5:
                        raise
                    sleep_s = min(2**attempts, 30)
                    log(f"  lock timeout on {spec.table}, retry in {sleep_s}s…")
                    time.sleep(sleep_s)

        if args.fake_django:
            fake_django_migration(conn, dry_run=args.dry_run)

        log("\nAll done.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
