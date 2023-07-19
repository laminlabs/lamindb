from typing import List, Union, overload  # noqa

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

        Delete files (delete the metadata record and the file in storage):

        >>> file = ln.select(File, id=file_id).one()
        >>> ln.delete(file)
        >>> # deleting the record occurs automatically
        >>> # you will be asked whether to delete the file in storage
        >>> # for more control, use:
        >>> file.delete(storage=True)

        Bulk delete via QuerySet:

        >>> ln.save(ln.Tag.from_values(["Tag1", "Tag2", "Tag3"], field="name"))
        >>> queryset = ln.Tag.select(name__icontains = "tag")
        >>> queryset.list()
        [Tag(id=o3FY3c5n, name=Tag2, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9), # noqa
        Tag(id=Qi3c4utq, name=Tag3, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9), # noqa
        Tag(id=CcFPLmpq, name=Tag1, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9)] # noqa
        >>> queryset.delete()
    """
    logger.warning("For efficient bulk delete, use `queryset.delete` instead")
    if isinstance(records, list):
        records = records
    elif isinstance(records, ORM):
        records = [records]
    for record in records:
        record.delete()
        logger.success(f"Deleted {colors.yellow(f'{record}')}")
