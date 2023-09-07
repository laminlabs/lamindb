from typing import List, Union, overload  # noqa

from lamin_utils import colors, logger
from lnschema_core import Registry


@overload
def delete(
    record: Registry,
) -> None:
    ...  # pragma: no cover


@overload
def delete(
    records: List[Registry],
) -> None:  # type: ignore
    ...  # pragma: no cover


def delete(  # type: ignore
    records: Union[Registry, List[Registry]],
) -> None:
    """Delete metadata records & files.

    Args:
        records: `Union[Registry, List[Registry]]` One or multiple records.

    Returns:
        `None`

    Examples:

        Delete by record:

        >>> experiment = ln.filter(Experiment, id=experiment_id).one()
        >>> ln.delete(experiment)

        Delete files (delete the metadata record and the file in storage):

        >>> file = ln.filter(File, id=file_id).one()
        >>> ln.delete(file)
        >>> # deleting the record occurs automatically
        >>> # you will be asked whether to delete the file in storage
        >>> # for more control, use:
        >>> file.delete(storage=True)

        Bulk delete via QuerySet:

        >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))
        >>> queryset = ln.ULabel.filter(name__icontains = "label")
        >>> queryset.list()
        [ULabel(id=o3FY3c5n, name=ULabel2, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9), # noqa
        ULabel(id=Qi3c4utq, name=ULabel3, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9), # noqa
        ULabel(id=CcFPLmpq, name=ULabel1, updated_at=2023-07-19 18:28:16, created_by_id=kmvZDIX9)] # noqa
        >>> queryset.delete()
    """
    logger.warning("for efficient bulk delete, use `queryset.delete` instead")
    if isinstance(records, list):
        records = records
    elif isinstance(records, Registry):
        records = [records]
    for record in records:
        record.delete()
        logger.success(f"deleted {colors.yellow(f'{record}')}")
