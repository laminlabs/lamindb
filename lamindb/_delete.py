from typing import List, Optional, Union, overload  # noqa

from django.db import transaction
from lamin_logger import colors, logger
from lnschema_core import ORM


@overload
def delete(
    record: ORM,
) -> None:
    ...


@overload
def delete(
    records: List[ORM],
) -> None:  # type: ignore
    ...


def delete(  # type: ignore
    records: Union[ORM, List[ORM]],
) -> None:
    """Delete metadata records & files.

    Args:
        records: `Union[ORM, List[ORM]]` One or multiple records.

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
    elif isinstance(records, ORM):
        records = [records]
    with transaction.atomic():
        for record in records:
            record.delete()
        logger.success(
            "Deleted"
            f" {colors.yellow(f'{len(records)}')} {ORM.__class__.__name__} records"
        )
