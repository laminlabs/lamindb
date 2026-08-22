from __future__ import annotations

from typing import TYPE_CHECKING

from django.db import ProgrammingError
from lamin_utils import logger
from lamindb_setup import settings as setup_settings

from ..core._settings import settings
from ..errors import NoWriteAccess
from .run import Run

WARNING_RUN_TRANSFORM = "no run & transform got linked, call `ln.track()` & re-run"
WARNING_NO_INPUT = "run input wasn't tracked, call `ln.track()` and re-run"

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .artifact import Artifact
    from .collection import Collection


def populate_recreating_run(dataset: Artifact | Collection, run: Run | None) -> None:
    if run is None:
        return
    if dataset.run is None:
        dataset.run = run
    elif dataset.run != run:
        # if this current run tracks this dataset already as
        # an input, then we don't track it as a recreating run to avoid a cycle
        # unfortunately we have to retrieve this information from the database
        # and can't cache it in dataset._input_of_run_id like we do for dataset._recreating_run_id
        # since recreating_runs is called in the object constructor
        # TODO: in the future we could rewrite the two network requests as a single using
        # plain SQL to optimize performance
        if dataset.input_of_runs.filter(id=run.id).exists():
            return
        dataset.recreating_runs.add(run)
        dataset._recreating_run_id = run.id


# also see current_run() in core._data
def get_run(run: Run | None) -> Run | None:
    from ..core._context import context
    from ..core._functions import get_current_tracked_run

    if run is None:
        run = get_current_tracked_run()
        if run is None:
            run = context.run
        if run is None and not settings.creation.artifact_silence_missing_run_warning:
            isettings = setup_settings.instance
            if not (isettings._is_clone or isettings.is_read_only_connection):
                logger.warning(WARNING_RUN_TRANSFORM)
    # suppress run by passing False
    elif not run:
        run = None
    return run


def is_valid_input(dataset: Artifact | Collection, run: Run) -> bool:
    is_valid = False
    # if a dataset is not yet saved it has dataset._state.db = None
    # then it can't be an input
    # we silently ignore because what will happen is that
    # the dataset either gets saved and then is tracked as an output
    # or it won't get saved at all
    if dataset._state.db == "default":
        # things are OK if the dataset is on the default db
        is_valid = True
    else:
        # dataset is on another db
        # we have to save the dataset into the current db with
        # the run being attached to a transfer transform
        logger.info(
            f"completing transfer to track {dataset.__class__.__name__}('{dataset.uid}') as input"
        )
        dataset.save()
        is_valid = True
    # avoid cycles: dataset can't be both input and output
    if dataset.run_id == run.id:
        logger.debug(
            f"not tracking {dataset} as input to run {run} because created by same run"
        )
        is_valid = False
    if run.id == getattr(dataset, "_recreating_run_id", None):
        logger.debug(
            f"not tracking {dataset} as input to run {run} because re-created in same run"
        )
        is_valid = False
    return is_valid


def track_run_inputs(
    dataset_or_datasets: Artifact
    | Iterable[Artifact]
    | Collection
    | Iterable[Collection],
    is_run_input: bool | Run | None = None,
    run: Run | None = None,
) -> None:
    """Links one or many datasets as inputs to a run.

    This function contains all validation logic to make decisions on whether a
    dataset qualifies as an input or not.
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
    # consider that dataset is an iterable of Data
    dataset_iter: Iterable[Artifact] | Iterable[Collection] = (
        [dataset_or_datasets]
        if isinstance(dataset_or_datasets, (Artifact, Collection))
        else dataset_or_datasets
    )
    input_datasets = []
    if run is not None:
        assert not run._state.adding, "Save the run before tracking its inputs."
        input_datasets = [
            dataset for dataset in dataset_iter if is_valid_input(dataset, run)
        ]
        input_datasets_ids = [dataset.id for dataset in input_datasets]
    if input_datasets:
        registry_str = input_datasets[0].__class__.__name__.lower()
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
        elif input_datasets:
            logger.debug(
                f"adding {registry_str} ids {input_datasets_ids} as inputs for run {run.id}"
            )
            track = True
    else:
        track = is_run_input
    if not track or not input_datasets:
        return None
    assert run is not None, "No run context set. Call `ln.track()`."
    if registry_str == "artifact":
        IsLink = run.input_artifacts.through
        links = [
            IsLink(run_id=run.id, artifact_id=dataset_id)
            for dataset_id in input_datasets_ids
        ]
    else:
        IsLink = run.input_collections.through
        links = [
            IsLink(run_id=run.id, collection_id=dataset_id)
            for dataset_id in input_datasets_ids
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
                dataset_space
                for dataset in input_datasets
                if (dataset_space := dataset.space) not in write_access_spaces
            }
            if (run_space := run.space) not in write_access_spaces:
                no_write_access_spaces.add(run_space)

            if not no_write_access_spaces:
                # if there are no unavailable spaces, then this should be due to locking
                locked_datasets = [
                    dataset
                    for dataset in input_datasets
                    if getattr(dataset, "is_locked", False)
                ]
                if run.is_locked:
                    locked_datasets.append(run)
                # if no unavailable spaces and no locked datasets, just raise the original error
                if not locked_datasets:
                    raise e
                no_write_msg = (
                    "It is not allowed to modify locked objects: "
                    + ", ".join(
                        r.__class__.__name__ + f"(uid={r.uid})" for r in locked_datasets
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
