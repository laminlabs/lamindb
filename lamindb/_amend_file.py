from pathlib import Path
from typing import Optional, Union, overload  # noqa

from lamindb_setup import settings as setup_settings
from lnschema_core.models import File, Run
from lnschema_core.types import DataLike, PathLike
from upath import UPath

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


def delete(file, storage: Optional[bool] = None) -> None:
    """Delete file, optionall from storage.

    Args:
        storage: `Optional[bool] = None` Indicate whether you want to delete the
        file in storage.
    """
    if storage is None:
        response = input(f"Are you sure you want to delete {file} from storage? (y/n)")
        if response == "y":
            delete_in_storage = True
    else:
        delete_in_storage = storage
    if delete_in_storage:
        filepath = file.path()
        delete_storage(filepath)
        logger.success(f"Deleted stored object {colors.yellow(f'{filepath}')}")
    file._delete_skip_storage()


def _delete_skip_storage(file, *args, **kwargs) -> None:
    super(File, file).delete(*args, **kwargs)


def save(file, *args, **kwargs) -> None:
    """Save the file to database & storage."""
    file._save_skip_storage(*args, **kwargs)
    from lamindb._save import check_and_attempt_clearing, check_and_attempt_upload

    exception = check_and_attempt_upload(file)
    if exception is not None:
        file._delete_skip_storage()
        raise RuntimeError(exception)
    exception = check_and_attempt_clearing(file)
    if exception is not None:
        raise RuntimeError(exception)


def _save_skip_storage(file, *args, **kwargs) -> None:
    if file.transform is not None:
        file.transform.save()
    if file.run is not None:
        file.run.save()
    super(File, file).save(*args, **kwargs)


def path(self) -> Union[Path, UPath]:
    """Path on storage."""
    from lamindb._file_access import filepath_from_file_or_folder

    return filepath_from_file_or_folder(self)


# likely needs an arg `key`
def replace(
    file,
    data: Union[PathLike, DataLike],
    run: Optional[Run] = None,
    format: Optional[str] = None,
) -> None:
    """Replace file content."""
    from lamindb._file import replace_file

    replace_file(file, data, run, format)


@overload
def __init__(
    file,
    data: Union[PathLike, DataLike],
    key: Optional[str] = None,
    name: Optional[str] = None,
    run: Optional[Run] = None,
):
    ...


@overload
def __init__(
    file,
    **kwargs,
):
    ...


def __init__(  # type: ignore
    file,
    *args,
    **kwargs,
):
    from lamindb._file import init_file

    init_file(file, *args, **kwargs)


File.backed = backed
File.stage = stage
File.load = load
File.delete = delete
File._delete_skip_storage = _delete_skip_storage
File.save = save
File._save_skip_storage = _save_skip_storage
File.replace = replace
File.__init__ = __init__
