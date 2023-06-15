from typing import List, Optional, Union, overload  # noqa

from lnschema_core import BaseORM

from ._logger import colors, logger


@overload
def delete(
    record: BaseORM,
) -> None:
    ...


@overload
def delete(
    records: List[BaseORM],
) -> None:  # type: ignore
    ...


def delete(  # type: ignore
    records: Union[BaseORM, List[BaseORM]],
) -> None:
    """Delete metadata records & files.

    Args:
        records: `Union[BaseORM, List[BaseORM]]` One or multiple records.

    Returns:
        `None`

    Examples:

        Delete by record:

        >>> experiment = ln.select(Experiment, id=experiment_id).one()
        >>> ln.delete(experiment)

        Delete files (delete the metadata record and the file in storage)

        >>> file = ln.select(File, id=file_id).one()
        >>> ln.delete(file)
        >>> # deleting the record occurs automatically
        >>> # you will be asked whether to delete the file in storage
        >>> # for more control, use:
        >>> file.delete(storage=True)

    """
    if isinstance(records, list):
        records = records
    elif isinstance(records, BaseORM):
        records = [records]
    for record in records:
        record.delete()
        logger.success(f"Deleted {colors.yellow(f'{record}')}")
