from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import ProgrammingError
from lamin_utils import logger
from lamindb_setup import settings as setup_settings

from ..core._settings import settings
from ..errors import NoWriteAccess
from .run import Run

WARNING_NO_INPUT = "run input wasn't tracked, call `ln.track()` and re-run"

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .artifact import Artifact
    from .collection import Collection


def is_valid_input(sqlrecord: Artifact | Collection, run: Run) -> bool:
    is_valid = False
    # if a sqlrecord is not yet saved it has sqlrecord._state.db = None
    # then it can't be an input
    # we silently ignore because what will happen is that
    # the sqlrecord either gets saved and then is tracked as an output
    # or it won't get saved at all
    if sqlrecord._state.db == "default":
        # things are OK if the sqlrecord is on the default db
        is_valid = True
    else:
        # sqlrecord is on another db
        # we have to save the sqlrecord into the current db with
        # the run being attached to a transfer transform
        logger.info(
            f"completing transfer to track {sqlrecord.__class__.__name__}('{sqlrecord.uid}') as input"
        )
        sqlrecord.save()
        is_valid = True
    # avoid cycles: sqlrecord can't be both input and output
    if sqlrecord.run_id == run.id:
        logger.debug(
            f"not tracking {sqlrecord} as input to run {run} because created by same run"
        )
        is_valid = False
    if run.id == getattr(sqlrecord, "_recreating_run_id", None):
        logger.debug(
            f"not tracking {sqlrecord} as input to run {run} because re-created in same run"
        )
        is_valid = False
    return is_valid


def track_run_inputs(
    sqlrecord_or_sqlrecords: (
        Artifact | Iterable[Artifact]
    ),  # can also be Collection | Iterable[Collection]
    is_run_input: bool | Run | None = None,
    run: Run | None = None,
) -> None:
    """Links one or many sqlrecords as inputs to a run.

    This function contains all validation logic to make decisions on whether a
    sqlrecord qualifies as an input or not.
    """
    if is_run_input is False:
        return None

    from ..core._context import context
    from ..core._functions import get_current_tracked_run
    from .artifact import Artifact
    from .collection import Collection

    if isinstance(is_run_input, Run):
        run = is_run_input
        is_run_input = True
    elif run is None:
        run = get_current_tracked_run()
        if run is None:
            run = context.run
    # consider that sqlrecord is an iterable of Data
    sqlrecord_iter: Iterable[Artifact] | Iterable[Collection] = (
        [sqlrecord_or_sqlrecords]
        if isinstance(sqlrecord_or_sqlrecords, (Artifact, Collection))
        else sqlrecord_or_sqlrecords
    )
    input_sqlrecords = []
    if run is not None:
        assert not run._state.adding, "Save the run before tracking its inputs."
        input_sqlrecords = [
            sqlrecord for sqlrecord in sqlrecord_iter if is_valid_input(sqlrecord, run)
        ]
        input_sqlrecords_ids = [sqlrecord.id for sqlrecord in input_sqlrecords]
    if input_sqlrecords:
        registry_str = input_sqlrecords[0].__class__.__name__.lower()
    # let us first look at the case in which the user does not
    # provide a boolean value for `is_run_input`
    # hence, we need to determine whether we actually want to
    # track a run or not
    track = False
    is_run_input = settings.track_run_inputs if is_run_input is None else is_run_input
    if is_run_input:
        if run is None:
            # Don't emit this warning when no global instance is configured.
            if setup_settings.is_configured:
                isettings = setup_settings.instance
                if not (isettings._is_clone or isettings.is_read_only_connection):
                    logger.warning(WARNING_NO_INPUT)
        elif input_sqlrecords:
            logger.debug(
                f"adding {registry_str} ids {input_sqlrecords_ids} as inputs for run {run.id}"
            )
            track = True
    else:
        track = is_run_input
    if not track or not input_sqlrecords:
        return None
    assert run is not None, "No run context set. Call `ln.track()`."
    if registry_str == "artifact":
        IsLink = run.input_artifacts.through
        links = [
            IsLink(run_id=run.id, artifact_id=sqlrecord_id)
            for sqlrecord_id in input_sqlrecords_ids
        ]
    else:
        IsLink = run.input_collections.through
        links = [
            IsLink(run_id=run.id, collection_id=sqlrecord_id)
            for sqlrecord_id in input_sqlrecords_ids
        ]
    try:
        IsLink.objects.bulk_create(links, ignore_conflicts=True)
    except ProgrammingError as e:
        if "new row violates row-level security policy" in str(e):
            instance = setup_settings.instance
            available_spaces = instance.available_spaces
            if available_spaces is None:
                raise NoWriteAccess(
                    f"You’re not allowed to write to the instance {instance.slug}.\n"
                    "Please contact administrators of the instance if you need write access."
                ) from None
            write_access_spaces = available_spaces["admin"] + available_spaces["write"]
            no_write_access_spaces = {
                sqlrecord_space
                for sqlrecord in input_sqlrecords
                if (sqlrecord_space := sqlrecord.space) not in write_access_spaces
            }
            if (run_space := run.space) not in write_access_spaces:
                no_write_access_spaces.add(run_space)

            if not no_write_access_spaces:
                # if there are no unavailable spaces, then this should be due to locking
                locked_sqlrecords = [
                    sqlrecord
                    for sqlrecord in input_sqlrecords
                    if getattr(sqlrecord, "is_locked", False)
                ]
                if run.is_locked:
                    locked_sqlrecords.append(run)
                # if no unavailable spaces and no locked sqlrecords, just raise the original error
                if not locked_sqlrecords:
                    raise e
                no_write_msg = (
                    "It is not allowed to modify locked sqlrecords: "
                    + ", ".join(
                        r.__class__.__name__ + f"(uid={r.uid})"
                        for r in locked_sqlrecords
                    )
                    + "."
                )
                raise NoWriteAccess(no_write_msg) from None

            if len(no_write_access_spaces) > 1:
                name_msg = ", ".join(
                    f"'{space.name}'" for space in no_write_access_spaces
                )
                space_msg = "spaces"
            else:
                name_msg = f"'{no_write_access_spaces.pop().name}'"
                space_msg = "space"
            raise NoWriteAccess(
                f"You’re not allowed to write to the {space_msg} {name_msg}.\n"
                f"Please contact administrators of the {space_msg} if you need write access."
            ) from None
        else:
            raise e


# backward compatibility
track_run_input = track_run_inputs
