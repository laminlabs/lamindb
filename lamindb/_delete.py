import traceback
from typing import List, Optional, Union, overload  # noqa

from lnschema_core import BaseORM, File, RunInput
from lnschema_core.models import storage_key_from_file

from lamindb.dev.storage import delete_storage

from ._logger import colors, logger
from ._select import select


@overload
def delete(
    record: BaseORM,
    delete_data_from_storage: Optional[bool] = None,
) -> None:
    ...


@overload
def delete(
    records: List[BaseORM],
    delete_data_from_storage: Optional[bool] = None,
) -> None:  # type: ignore
    ...


@overload
def delete(
    entity: BaseORM,
    delete_data_from_storage: Optional[bool] = None,
    **fields,
) -> None:  # type: ignore
    ...


def delete(  # type: ignore
    record: Union[BaseORM, List[BaseORM]],
    delete_data_from_storage: Optional[bool] = None,
    **fields,
) -> None:
    """Delete metadata records & files.

    Args:
        record: One or multiple records as instances of `SQLModel`.
        delete_data_from_storage: Whether to delete data from storage.

    Returns:
        `None`

    Examples:

        Delete by record:

        >>> experiment = ln.select(Experiment, id=experiment_id).one()
        >>> ln.delete(experiment)

        Delete by fields:

        >>> ln.delete(Experiment, id=experiment_id)
        >>> # the result of is equivalent to 1)

        Delete files (delete the metadata record and the file in storage)

        >>> file = ln.select(File, id=file_id).one()
        >>> # deleting the metadata record occurs automatically
        >>> # you will be asked whether to delete the file from storage
        >>> # or pass boolean values to `delete_data_from_storage`
        >>> ln.delete(file, delete_data_from_storage)

    """
    if isinstance(record, list):
        records = record
    elif isinstance(record, BaseORM):
        records = [record]
    else:
        model = record
        results = select(model, **fields).one_or_none()
        if results is None:
            return results
        else:
            records = [results]

    for record in records:
        storage_key = None
        if isinstance(record, File):
            # save storage key before deleting the record
            # after the deletion file.id is None
            storage_key = storage_key_from_file(record)
            # delete run_ins related to the file that's to be deleted
            run_ins = select(RunInput, file_id=record.id).all()
            for run_in in run_ins:
                run_in.delete()
        try:
            record.delete()
            logger.success(
                f"Deleted {colors.yellow(f'row {record}')} in"
                f" {colors.blue(f'table {type(record).__name__}')}."
            )
        except Exception:
            traceback.print_exc()
        if isinstance(record, File):
            if delete_data_from_storage is None:
                # ask to confirm deleting data from storage
                delete_dialog = (
                    "Confirm Delete: Are you sure you want to delete"
                    f" object {storage_key} from storage? (y/n)"
                )
                decide = input(f"   {delete_dialog}")
            else:
                decide = "y" if delete_data_from_storage else "n"

            if decide not in ("y", "Y", "yes", "Yes", "YES"):
                continue
            try:
                delete_storage(storage_key)  # type: ignore
                logger.success(
                    f"Deleted {colors.yellow(f'object {storage_key}')} from storage."
                )
            except Exception:
                traceback.print_exc()
