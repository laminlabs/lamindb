# Schema hash lookup implementation plan

## Decision

Keep schema hashes structural and make constructor lookup respect an explicitly
supplied schema type. Identical schemas may have the same hash under different
types without moving an existing record into the requested type.

References: [issue 1285](https://github.com/pfizer-collab/laminlabs/issues/1285),
[Slack discussion](https://laminlabs.slack.com/archives/C07D9G0BN1J/p1794491657105649),
and [PR 4078](https://github.com/laminlabs/lamindb/pull/4078).

The baseline test on main `a1752c4b` fails after creating the second schema:
reloading the first shows `type_id=2` instead of `type_id=1`. The constructor finds
the original record by hash, copies its identity, and applies the second type.

## Lookup contract

| Constructor argument     | Search scope                                    | Multiple matches             |
| ------------------------ | ----------------------------------------------- | ---------------------------- |
| Omitted `type` (`UNSET`) | All types, including root                       | Raise `MultipleResultsFound` |
| Explicit `type=None`     | Root schemas only                               | Raise `MultipleResultsFound` |
| Explicit saved type      | That immediate type only, excluding descendants | Raise `MultipleResultsFound` |

Continue excluding trashed schemas. Keep `.one_or_none()`; never select an
arbitrary match. No match proceeds through the existing type-aware name lookup
and then normal creation. `lookup=False` continues to bypass only hash lookup.
Type records (`is_type=True`) continue to have no structural hash.

## Where to change the code

1. **`lamindb/models/schema.py`, `Schema.__init__`:** build the existing hash
   queryset, then add `type=type` only when the argument is not `UNSET`. This is
   the smallest natural change and preserves existing hash calculation, branch
   exclusion, database selection, and name fallback behavior. Document the three
   argument states in the constructor documentation.
2. **Built-in UID assignment in the same constructor:** reserve a known schema's
   canonical UID for root construction. Explicitly typed copies receive generated
   UIDs. If a canonical UID is already occupied (for example, its original record
   was moved to a type or trashed), a new root record also receives a generated
   UID. Check UID occupancy including trash because the database uniqueness
   constraint includes it. Do not rewrite existing UIDs or hashes.
3. **`Schema.save()`:** when structural changes would duplicate a hash, check
   only the record's immediate type, including root. Keep the existing rejection
   of structural edits that duplicate a schema in that same scope. Improve the
   error wording so it describes the actual hash conflict rather than only feature
   order.

These are schema-specific rules; do not change shared name lookup or introduce a
new general lookup abstraction for two short queries. No database migration or
rehashing is needed.

## Affected methods and boundaries

- **`Schema.from_values()` and `Schema.from_dataframe()`:** call the constructor
  without a type. They still reuse a unique match across types and now raise if
  the newly supported duplicates make that lookup ambiguous.
- **`curators/core.py`, `annotate_artifact()`:** inferred feature sets still omit
  `type`, as agreed in Slack. No automatic inheritance from `artifact.schema.type`
  is introduced. Duplicate matches raise; an explicit schema does not yet define
  type propagation for every inferred feature set.
- **`FeatureManager` validation and `_add_from()`:** inferred/transfer feature-set
  lookup remains unscoped. The direct hash lookup already uses `.one_or_none()`,
  and the other path calls `Schema.from_values()`. Both reject ambiguity rather
  than choosing a type. Transfer redesign and the pre-existing UID-copy behavior
  are outside this patch; no new `.first()` fallback is added.
- **Composite schemas:** slot hashes remain based on slot names and child hashes.
  Equal-hash children under different types can produce equal-hash parents.
  Different explicit parent types separate the parent records. Within one parent
  type, hash reuse keeps the existing parent's child links, even if a different
  equal-hash child was supplied. This accepted structural-equivalence tradeoff is
  documented and covered; preserving exact child identity would require a
  separate design.
- **Type-only moves:** preserve the hash and UID. No new database uniqueness or
  automatic merge policy is introduced. A move can create multiple records with
  the same `(type, hash)`; subsequent constructor lookups raise. This is distinct
  from the existing guard on structural edits and is tested explicitly.
- **Artifacts, features, records, and collections:** their hashes and constructor
  matching rules remain unchanged.

## Edge cases and verification plan

| Case                                                          | Expected behavior                                      |
| ------------------------------------------------------------- | ------------------------------------------------------ |
| Original issue, same name and structure under different types | Distinct UIDs; first type preserved; equal hashes      |
| Same structure under the same explicit type                   | Reuse existing record                                  |
| Omitted type with one typed match                             | Reuse without changing its type                        |
| Omitted type with multiple matches, even with a matching name | Raise before updating any record                       |
| Explicit `None` with only typed matches                       | Create a root record; never adopt a typed record       |
| Root plus typed duplicates                                    | Explicit scope resolves each correctly                 |
| Parent and child schema types                                 | Match only the immediate type                          |
| Trashed match                                                 | Excluded from lookup                                   |
| Unsaved type                                                  | Reject rather than reuse and mutate a saved record     |
| `lookup=False`                                                | No hash lookup or ambiguity error                      |
| Built-in schema created root-first or typed-first             | Root canonical UID when available; distinct typed UIDs |
| Canonical UID occupied by a moved or trashed record           | New root uses a distinct generated UID                 |
| Structural edit matches another type                          | Allow; both records remain distinct                    |
| Structural edit matches the same type, including root         | Reject before persistence                              |
| Type-only move, including into an occupied scope              | Hash and UID unchanged; ambiguous future lookup raises |
| Composite parent reuse                                        | Preserve structural hashes and existing child links    |

Use the existing isolated SQLite pytest fixture and small synthetic features.
Keep cleanup in `finally`. Run the original regression, the added edge cases,
and relevant existing schema tests. Run broader schema/annotation tests where
their local data dependencies are available; distinguish environment failures
from regressions. No playground or production data is needed for these unit tests.

## Commit and delivery

Commit the implementation, tests, and this document on
`fix/schema-hash-type-lookup`. Use `[skip ci]` in the commit message and verify no
new Actions run is associated with that commit. The current workflows use
`pull_request` for branch updates, which supports GitHub's
[skip instruction](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/skip-workflow-runs).
Update the existing draft PR with the implementation and local verification.

## Verification results

Implemented the three changes above in `lamindb/models/schema.py`. The original
regression test remains unchanged. Additional coverage is in
`tests/pydata/test_schema_hash_lookup.py`; its fixture rolls back all test records.

### Before the fix

The original regression failed locally on SQLite (`2 == 1`) and in both PR CI
jobs at commit `f2c96fb`:

- [SQLite baseline](https://github.com/laminlabs/lamindb/actions/runs/43267120905/job/131353674724): `312 == 311`.
- [PostgreSQL baseline](https://github.com/laminlabs/lamindb/actions/runs/43267120905/job/131353674524): `460 == 459`.

All failures are at the assertion that the first schema retains its original
`type_id` after saving the second schema.

### After the fix

- Original regression plus the initial edge-case suite: **15 passed**.
- Existing schema tests, the expanded edge-case suite, and dataframe annotation
  tests together: **96 passed, 6 skipped** in 38.63 seconds. Skips are optional
  dependency cases in the annotation tests.
- Tested locally with Python 3.12.3, SQLite, this worktree, and the dependency
  submodule revisions pinned by main `a1752c4b`. Settings and caches use dedicated
  temporary directories; the suite removes its temporary database and storage.
- The fix has not been rerun on PostgreSQL. CI is intentionally skipped for the
  implementation push as requested. Cross-instance transfer paths were reviewed
  but were not changed or executed in this verification.

The broader test invocation, from the worktree:

```sh
env \
  LAMIN_SETTINGS_DIR=/private/tmp/lamindb-schema-hash-broader-settings \
  LAMIN_CACHE_DIR=/private/tmp/lamindb-schema-hash-repro-cache \
  BIONTY_ROOT_DIR=/private/tmp/lamindb-schema-hash-repro-bionty \
  PYTHONPATH="$PWD:$PWD/sub/lamindb-setup:$PWD/sub/bionty:$PWD/sub/lamin-cli" \
  LAMINDB_TEST_DB_VENDOR=sqlite \
  /Users/joelottosson/Documents/programming/laminlabs/lamindb/.venv/bin/python \
  -m pytest tests/pydata/test_schema.py tests/pydata/test_schema_hash_lookup.py \
  tests/pydata/test_artifact_dataframe_with_curation.py -q --tb=short
```

The isolated Bionty root contains a copy of `sub/bionty/bionty/base/sources.yaml`.
No runtime dependencies were changed. Repository commit checks cover formatting,
lint, and type checking before the final commit.
