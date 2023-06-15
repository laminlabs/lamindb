from pathlib import Path
from typing import Optional

from lamindb_setup import settings as setup_settings
from lnschema_core import File
from lnschema_core.types import DataLike

from lamindb._context import context
from lamindb._file_access import filepath_from_file_or_folder
from lamindb.dev.storage import delete_storage, load_to_memory
from lamindb.dev.storage.object._anndata_accessor import AnnDataAccessor

from ._logger import colors, logger
from ._settings import settings

File.__doc__ = """Files: data artifacts.

Args:
   data: `Union[PathLike, DataLike]` A file path or an in-memory data
      object (`DataFrame`, `AnnData`) to serialize. Can be a cloud path, e.g.,
      `"s3://my-bucket/my_samples/my_file.fcs"`.
   key: `Optional[str] = None` A storage key: a relative filepath within the
      current default storage, e.g., `"my_samples/my_file.fcs"`.
   name: `Optional[str] = None` A name. Defaults to the filename.
   run: `Optional[Run] = None` The run that created the file.

Track where files come from by passing the generating :class:`~lamindb.Run`.

Often, files store jointly measured observations of features: track them
with :class:`~lamindb.FeatureSet`.

If files have corresponding representations in storage and memory, LaminDB
makes some configurable default choices (e.g., serialize a `DataFrame` as a
`.parquet` file).

.. admonition:: Examples for storage-memory correspondence

   Listed are typical `suffix` values & in memory data objects.

   - Table: `.csv`, `.tsv`, `.parquet`, `.ipc`
     ⟷ `pd.DataFrame`, `polars.DataFrame`
   - Annotated matrix: `.h5ad`, `.h5mu`, `.zrad` ⟷ `AnnData`, `MuData`
   - Image: `.jpg`, `.png` ⟷ `np.ndarray`, ...
   - Array: zarr directory, TileDB store ⟷ zarr loader, TileDB loader
   - Fastq: `.fastq` ⟷ /
   - VCF: `.vcf` ⟷ /
   - QC: `.html` ⟷ /

"""


def backed(file: File, is_run_input: Optional[bool] = None) -> AnnDataAccessor:
    """Return a cloud-backed data object to stream."""
    _track_run_input(file, is_run_input)
    if file.suffix not in (".h5ad", ".zrad", ".zarr"):
        raise ValueError("File should have an AnnData object as the underlying data")
    return AnnDataAccessor(file)


def _track_run_input(file: File, is_run_input: Optional[bool] = None):
    if is_run_input is None:
        if context.run is not None:
            logger.hint("Track this file as a run input by passing `is_run_input=True`")
        track_run_input = settings.track_run_inputs_upon_load
    else:
        track_run_input = is_run_input
    if track_run_input:
        if context.run is None:
            raise ValueError(
                "No global run context set. Call ln.context.track() or link input to a"
                " run object via `run.inputs.append(file)`"
            )
        if not file.input_of.contains(context.run):
            context.run.save()
            file.input_of.add(context.run)


def load(
    file: File, is_run_input: Optional[bool] = None, stream: bool = False
) -> DataLike:
    """Stage and load to memory.

    Returns in-memory representation if possible, e.g., an `AnnData` object
    for an `h5ad` file.
    """
    if stream and file.suffix not in (".h5ad", ".zrad", ".zarr"):
        raise ValueError(
            "For streaming, file should have an AnnData object as the underlying data"
        )
    _track_run_input(file, is_run_input)
    return load_to_memory(filepath_from_file_or_folder(file), stream=stream)


def stage(file: File, is_run_input: Optional[bool] = None) -> Path:
    """Update cache from cloud storage if outdated.

    Returns a path to a locally cached on-disk object (say, a
    `.jpg` file).
    """
    if file.suffix in (".zrad", ".zarr"):
        raise RuntimeError("zarr object can't be staged, please use load() or stream()")
    _track_run_input(file, is_run_input)
    return setup_settings.instance.storage.cloud_to_local(
        filepath_from_file_or_folder(file)
    )


def delete(file, *args, **kwargs) -> None:
    delete_in_storage = False
    if "storage" in kwargs:
        delete_in_storage = kwargs.pop("storage")
    else:
        response = input(
            f"Are you sure you want to delete file {file} from storage? (y/n)"
        )
        if response == "y":
            delete_in_storage = True
    file._delete_skip_storage(*args, **kwargs)
    if delete_in_storage:
        filepath = file.path()
        delete_storage(filepath)
        logger.success(f"Deleted stored file at {colors.yellow(f'{filepath}')}")


def _delete_skip_storage(file, *args, **kwargs) -> None:
    super(File, file).delete(*args, **kwargs)


File.backed = backed
File.stage = stage
File.load = load
File.delete = delete
File._delete_skip_storage = _delete_skip_storage
