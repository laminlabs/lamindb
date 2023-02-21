from pathlib import Path
from typing import Optional, Union

from lndb import settings as setup_settings
from lndb.dev import UPath
from lnschema_core import DObject as lns_DObject
from lnschema_core import Run, Storage
from sqlalchemy.exc import NoResultFound

from .dev.db._select import select


def get_dfolder_kwargs_from_data(
    folder: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    source: Optional[Run] = None,
):
    folderpath = UPath(folder)
    storage_root = str(setup_settings.instance.storage.root)
    # backwards compatibility for instances created with cloudpathlib
    try:
        storage = select(Storage, root=storage_root).one()  # noqa
    except NoResultFound:
        storage = select(Storage, root=storage_root[:-1]).one()  # noqa

    dobjects = []
    for f in folderpath.glob("**/*"):
        if f.is_file():
            dobjects.append(lns_DObject(f, source=source))

    dfolder_kwargs = dict(
        name=folderpath.name if name is None else name,
        dobjects=dobjects,
    )
    return dfolder_kwargs
