from typing import List, Optional, Union, overload  # noqa

from lnschema_core import BaseORM

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
        record.delete()
        logger.success(
            f"Deleted {colors.yellow(f'row {record}')} in"
            f" {colors.blue(f'table {type(record).__name__}')}"
        )
