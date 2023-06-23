from itertools import islice
from pathlib import Path
from typing import Optional, Union, overload  # noqa

from lamin_logger import colors, logger
from lamindb_setup import settings as setup_settings
from lnschema_core.models import File, Run
from lnschema_core.types import DataLike, PathLike
from upath import UPath

from lamindb._context import context
from lamindb._file import from_dir, init_file, replace_file
from lamindb._file_access import filepath_from_file
from lamindb.dev._settings import settings
from lamindb.dev.storage import delete_storage, load_to_memory

try:
    from lamindb.dev.storage._backed_access import AnnDataAccessor, BackedAccessor
except ImportError:

    class AnnDataAccessor:  # type: ignore
        pass

    class BackedAccessor:  # type: ignore
        pass


File.__doc__ = """Files: data artifacts.

Args:
   data: `Union[PathLike, DataLike]` A file path or an in-memory data
      object (`DataFrame`, `AnnData`) to serialize. Can be a cloud path, e.g.,
      `"s3://my-bucket/my_samples/my_file.fcs"`.
   key: `Optional[str] = None` A storage key: a relative filepath within the
      current default storage, e.g., `"my_samples/my_file.fcs"`.
   name: `Optional[str] = None` A name or title. Useful if key is auto-generated.
   run: `Optional[Run] = None` The run that created the file, gets auto-linked
       if `ln.track()` was called.

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

.. note::

    In some cases (`.zarr`), a `File` is present as many small objects in what
    appears to be a "folder" in storage. Hence, we often refer to files as data
    artifacts.

"""


def backed(
    file: File, is_run_input: Optional[bool] = None
) -> Union[AnnDataAccessor, BackedAccessor]:
    """Return a cloud-backed data object to stream."""
    suffixes = (".h5", ".hdf5", ".h5ad", ".zrad", ".zarr")
    if file.suffix not in suffixes:
        raise ValueError(
            "File should have a zarr or h5 object as the underlying data, please use"
            " one of the following suffixes for the object name:"
            f" {', '.join(suffixes)}."
        )
    _track_run_input(file, is_run_input)
    from lamindb.dev.storage._backed_access import backed_access

    return backed_access(file)


def _track_run_input(file: File, is_run_input: Optional[bool] = None):
    if is_run_input is None:
        if context.run is not None and not settings.track_run_inputs:
            logger.hint("Track this file as a run input by passing `is_run_input=True`")
        track_run_input = settings.track_run_inputs
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
    _track_run_input(file, is_run_input)
    return load_to_memory(filepath_from_file(file), stream=stream)


def stage(file: File, is_run_input: Optional[bool] = None) -> Path:
    """Update cache from cloud storage if outdated.

    Returns a path to a locally cached on-disk object (say, a
    `.jpg` file).
    """
    if file.suffix in (".zrad", ".zarr"):
        raise RuntimeError("zarr object can't be staged, please use load() or stream()")
    _track_run_input(file, is_run_input)
    return setup_settings.instance.storage.cloud_to_local(filepath_from_file(file))


def delete(file, storage: Optional[bool] = None) -> None:
    """Delete file, optionall from storage.

    Args:
        storage: `Optional[bool] = None` Indicate whether you want to delete the
        file in storage.

    Example:

    For any `File` object `file`, call:

    >>> file.delete(storage=True)  # storage=True auto-confirms deletion in storage
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
    from lamindb._file_access import filepath_from_file

    return filepath_from_file(self)


# adapted from: https://stackoverflow.com/questions/9727673/list-directory-tree-structure-in-python  # noqa
@classmethod  # type: ignore
def tree(
    cls: File,
    prefix: Optional[str] = None,
    *,
    level: int = -1,
    limit_to_directories: bool = False,
    length_limit: int = 1000,
):
    """Given a prefix, print a visual tree structure of files."""
    space = "    "
    branch = "│   "
    tee = "├── "
    last = "└── "

    if prefix is None:
        dir_path = settings.storage
    else:
        dir_path = settings.storage / prefix
    files = 0
    directories = 0

    def inner(dir_path: Union[Path, UPath], prefix: str = "", level=-1):
        nonlocal files, directories
        if not level:
            return  # 0, stop iterating
        stripped_dir_path = dir_path.as_posix().rstrip("/")
        # do not iterate through zarr directories
        if stripped_dir_path.endswith((".zarr", ".zrad")):
            return
        # this is needed so that the passed folder is not listed
        contents = [
            i
            for i in dir_path.iterdir()
            if i.as_posix().rstrip("/") != stripped_dir_path
        ]
        if limit_to_directories:
            contents = [d for d in contents if d.is_dir()]
        pointers = [tee] * (len(contents) - 1) + [last]
        for pointer, path in zip(pointers, contents):
            if path.is_dir():
                yield prefix + pointer + path.name
                directories += 1
                extension = branch if pointer == tee else space
                yield from inner(path, prefix=prefix + extension, level=level - 1)
            elif not limit_to_directories:
                yield prefix + pointer + path.name
                files += 1

    folder_tree = f"{dir_path.name}"
    iterator = inner(dir_path, level=level)
    for line in islice(iterator, length_limit):
        folder_tree += f"\n{line}"
    if next(iterator, None):
        folder_tree += f"... length_limit, {length_limit}, reached, counted:"
    print(folder_tree)
    print(f"\n{directories} directories" + (f", {files} files" if files else ""))


# likely needs an arg `key`
def replace(
    file,
    data: Union[PathLike, DataLike],
    run: Optional[Run] = None,
    format: Optional[str] = None,
) -> None:
    """Replace file content.

    Args:
        data: `Union[PathLike, DataLike]` A file path or an in-memory data
            object (`DataFrame`, `AnnData`).
        run: `Optional[Run] = None` The run that created the file, gets
            auto-linked if `ln.track()` was called.

    Examples:

    Say we made a change to the content of a file (e.g., edited the image
    `paradisi05_laminopathic_nuclei.jpg`).

    This is how we replace the old file in storage with the new file:

    >>> file.replace("paradisi05_laminopathic_nuclei.jpg")
    >>> file.save()

    Note that this neither changes the storage key nor the filename.

    However, it will update the suffix if the file type changes.
    """
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
File.path = path
File.from_dir = from_dir
File.tree = tree
