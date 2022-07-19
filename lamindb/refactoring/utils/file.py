from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath


def create_dir_if_not_exists(path: Union[Path, CloudPath]):
    if not path.exists():
        path.mkdir(parents=True)
