# Tested in lamin-cli (tests/core/test_create_switch_delete_list_settings.py::test_merge*).
from __future__ import annotations

from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from django.apps import apps
from django.db import connection
from django.db.utils import DatabaseError
from lamin_utils import logger

if TYPE_CHECKING:
    from lamindb.models import Branch


def merge(branch: str | Branch) -> None:
    """Merge a branch into the current branch.

    All `SQLRecord` objects that have `branch_id` equal to the source branch's id
    are updated to the current branch's id.

    Find more info in the :class:`~lamindb.Branch` document.

    Args:
        branch: The source branch to merge from. Accepts a `name`, a `uid`, or the `Branch` object.

    Raises:
        DoesNotExist: If the branch does not exist.
    """
    from lamindb import Branch, Q
    from lamindb.errors import ObjectDoesNotExist

    from ..models import SQLRecord
    from ..models._is_versioned import IsVersioned, reconcile_is_latest_within_branch

    if isinstance(branch, Branch):
        source = branch
        if source._state.adding:
            raise ObjectDoesNotExist("Branch must be saved.")
    else:
        source = Branch.filter(Q(name=branch) | Q(uid=branch)).one_or_none()
        if source is None:
            raise ObjectDoesNotExist(f"Branch '{branch}' not found.")

    current = ln_setup.settings.branch
    if current.id == source.id:
        logger.important("already on branch, nothing to merge")
        return

    models = [
        m
        for m in apps.get_models()
        if issubclass(m, SQLRecord) and not m._meta.abstract
    ]
    if not models:
        return

    vendor = connection.vendor
    quoted_tables = [connection.ops.quote_name(m._meta.db_table) for m in models]

    with connection.cursor() as cursor:
        if vendor == "postgresql":
            # Single round-trip: one multi-statement execute
            statements = [
                f"UPDATE {tbl} SET branch_id = %s WHERE branch_id = %s"
                for tbl in quoted_tables
            ]
            sql = "BEGIN; " + "; ".join(statements) + "; COMMIT;"
            params = [current.id, source.id] * len(quoted_tables)
            try:
                cursor.execute(sql, params)
            except DatabaseError as e:
                logger.error(f"Merge failed: {e}")
                raise
        else:
            # SQLite: execute() runs only the first statement; run each UPDATE
            # in a loop (same connection, so still one transaction if we're inside
            # a transaction or use autocommit-off).
            from django.db import transaction

            with transaction.atomic():
                for tbl in quoted_tables:
                    # Django uses %s; SQLite backend converts to ?
                    cursor.execute(
                        f"UPDATE {tbl} SET branch_id = %s WHERE branch_id = %s",
                        [current.id, source.id],
                    )

    versioned_models = [m for m in models if issubclass(m, IsVersioned)]
    for model in versioned_models:
        reconcile_is_latest_within_branch(model, branch_id=current.id)

    source._status_code = -1  # merged
    source.save(update_fields=["_status_code"])
    logger.important(f"merged branch '{source.name}' into '{current.name}'")
