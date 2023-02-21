from pathlib import Path
from typing import Optional, Union

from lndb import settings as setup_settings
from lndb.dev import UPath
from lnschema_core import DObject as lns_DObject
from lnschema_core import Run, Storage
from sqlalchemy.exc import NoResultFound

from ._record import get_run
from .dev.db._select import select


def get_dfolder_kwargs_from_data(
    folder: Union[Path, UPath, str],
    *,
    name: Optional[str] = None,
    source: Optional[Run] = None,
):
    run = get_run(source)
    folderpath = UPath(folder)
    storage_root = str(setup_settings.instance.storage.root)
    # backwards compatibility for instances created with cloudpathlib
    try:
        storage = select(Storage, root=storage_root).one()
    except NoResultFound:
        storage = select(Storage, root=storage_root[:-1]).one()

    dobjects = []
    for f in folderpath.glob("**/*"):
        if f.is_file():
            dobjects.append(lns_DObject(f))

    dfolder_kwargs = dict(
        name=name,
        run_id=run.id,
        storage_id=storage.id,
        source=run,
        dobjects=dobjects,
    )
    return dfolder_kwargs
