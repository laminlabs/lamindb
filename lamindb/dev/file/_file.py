import shutil
from pathlib import Path
from typing import Union

import anndata as ad
import pandas as pd
from cloudpathlib import CloudPath
from lndb_setup import settings

READER_FUNCS = {".csv": pd.read_csv, ".h5ad": ad.read, ".feather": pd.read_feather}


def store_file(filepath: Union[str, Path], filekey: str):
    """Store arbitrary file."""
    storage_path = settings.instance.storage.key_to_filepath(filekey)
    if isinstance(storage_path, CloudPath):
        storage_path.upload_from(filepath)
    else:
        shutil.copyfile(filepath, storage_path)


def load_to_memory(filepath: Union[str, Path]):
    filepath = Path(filepath)
    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        raise NotImplementedError
    else:
        return reader(filepath)
