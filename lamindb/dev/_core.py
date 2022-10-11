from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from lndb_setup import settings
from lnschema_core import dobject


def get_name_suffix_from_filepath(filepath: Union[Path, CloudPath]):
    suffix = "".join(filepath.suffixes)
    name = filepath.name.replace(suffix, "")
    return name, suffix


def storage_key_from_dobject(dobj: dobject):
    return f"{dobj.id}{dobj.suffix}"


def filepath_from_dobject(dobj: dobject):
    storage_key = storage_key_from_dobject(dobj)
    filepath = settings.instance.storage.key_to_filepath(storage_key)
    return filepath
