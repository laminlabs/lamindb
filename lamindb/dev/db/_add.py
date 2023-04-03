from functools import partial
from pathlib import Path, PurePath
from typing import List, Optional, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lamin_logger import logger
from lndb import settings as setup_settings
from lndb_storage import UPath, delete_storage, store_object, write_adata_zarr
from lndb_storage._file import print_hook
from lnschema_core import File, Folder
from pydantic.fields import ModelPrivateAttr
from sqlalchemy.orm.attributes import set_attribute

from .._docs import doc_args
from ._core import file_to_sqm, get_session_from_kwargs
from ._select import select

add_docs = """
Insert or update data records.

Inserts a new :term:`record` if the corresponding row doesn't exist.
Updates the corresponding row with the record if it exists.

To update a row, query it with `.get` or `.select` and modify it before
passing it to `add`.

Guide: :doc:`/guide/add-delete`.

Examples:

1) Add a record (by passing a record)
>>> ln.add(lns.Pipeline(name="My pipeline", v="1"))
Pipeline(id="0Cb86EZj", name="My pipeline", v="1", ...)

2) Update an existing record
>>> pipeline = ln.select(lns.Pipeline, id="0Cb86EZj").one()
>>> pipeline.name = "New name"
>>> ln.add(experiment)
Pipeline(id="0Cb86EZj", name="New name", v="1", ...)

3) Add a record if not exist in the DB
>>> # add a record if the metadata combination is not already exist in the DB
>>> # if exists, returns the existing record from the DB
>>> ln.add(lns.Pipeline, name="My pipeline", v="1")
>>> # under the hood, this runs a query first based on passed fields
>>> # equivalent to the following:
>>> pipeline = ln.select(lns.Pipeline, name="My pipeline", v="1").one_or_none()
>>> if pipeline is None:
>>>     ln.add(pipeline)

Args:
    record: One or multiple records as instances of `SQLModel`.
"""


@overload
def add(record: sqm.SQLModel) -> sqm.SQLModel:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def add(records: List[sqm.SQLModel]) -> List[sqm.SQLModel]:  # type: ignore
    ...


@overload
def add(  # type: ignore
    entity: sqm.SQLModel, **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    ...


@doc_args(add_docs)
def add(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]], **fields
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """{}"""  # noqa
    session = get_session_from_kwargs(fields)
    if isinstance(record, list):
        records = record
    elif isinstance(record, sqm.SQLModel):
        records = [record]
    else:
        model = file_to_sqm(record)
        results = select(model, **fields).one_or_none()
        if results is None:
            records = [model(**fields)]
        else:
            logger.info(
                f"An existing {results.__class__.__name__} record is found in the DB:"
            )
            return results

    if session is None:  # assume global session
        session = setup_settings.instance.session()
        setup_settings.instance._cloud_sqlite_locker.lock()
        close = True
    else:
        close = False

    # commit metadata to database
    db_error = None
    for record in records:
        write_objectkey(record)
        # the following ensures that queried objects (within __init__)
        # behave like queried objects, only example right now: Run
        if hasattr(record, "_ln_identity_key") and record._ln_identity_key is not None:
            record._sa_instance_state.key = record._ln_identity_key
        session.add(record)
    try:
        session.commit()
    except Exception as e:
        db_error = e

    # upload data objects to storage
    added_records = []
    if db_error is None:
        added_records, upload_error = upload_committed_records(records, session)

    if close:
        session.close()
        setup_settings.instance._update_cloud_sqlite_file()
        setup_settings.instance._cloud_sqlite_locker.unlock()

    error = db_error or upload_error
    if error is not None:
        error_message = prepare_error_message(records, added_records, error)
        raise RuntimeError(error_message)
    elif len(added_records) > 1:
        return added_records
    else:
        return added_records[0]


def upload_committed_records(records, session):
    """Upload records in a list of database-commited records to storage.

    If any upload fails, subsequent records are cleaned up from the DB.
    """
    # make sure ALL records are up-to-date to enable accurate comparison
    # during metadata cleanup
    error = None
    added_records = []
    for record in records:
        session.refresh(record)

    # upload data objects
    for record in records:
        if isinstance(record, File) and hasattr(record, "_local_filepath"):
            try:
                upload_data_object(record)
            except Exception as e:
                error = e
                break
        added_records += [record]

    # clear old files on update
    for record in added_records:
        if isinstance(record, File) and hasattr(record, "_clear_storagekey"):
            try:
                if record._clear_storagekey is not None:
                    delete_storage(record._clear_storagekey)
                    record._clear_storagekey = None
            except Exception as e:
                error = e

    # clean up metadata for objects not uploaded to storage
    if error is not None:
        for record in records:
            if record not in added_records:
                session.delete(record)
        session.commit()

    return added_records, error


def prepare_error_message(records, added_records, error) -> str:
    if len(records) == 1 or len(added_records) == 0:
        error_message = (
            "An error occured. No entries were uploaded or committed"
            " to the database. See error message below.\n\n"
        )
    else:
        error_message = (
            "An error occured. The following entries have been"
            " successfully uploaded and committed to the database:\n"
        )
        for record in added_records:
            error_message += (
                f"- {', '.join(record.__repr__().split(', ')[:3]) + ', ...)'}\n"
            )
        error_message += "\nSee error message below.\n\n"
    error_message += f"{str(error)}"
    return error_message


def local_instance_storage_matches_local_parent(file: File):
    storage = setup_settings.instance.storage

    if file._local_filepath is not None:
        path = file._local_filepath
    else:
        path = file.path()

    if not isinstance(path, UPath):
        path = path.resolve()

    parents = [str(p) for p in path.parents]
    return str(storage.root) in parents


def get_storage_root_and_root_str(
    root: Optional[Union[Path, UPath]] = None
) -> Tuple[Union[Path, UPath], str]:
    if root is None:
        root = setup_settings.instance.storage.root
    root_str = root.as_posix()
    if isinstance(root, UPath):
        root_str = root_str.rstrip("/")
    return root, root_str


def filepath_to_relpath(
    root: Union[PurePath, Path], root_str: str, filepath: Union[Path, UPath]
) -> Union[PurePath, Path]:
    """Filepath to relative path of the root."""
    if isinstance(root, UPath):
        relpath = PurePath(filepath.as_posix().replace(root_str, ""))
    else:
        relpath = filepath.resolve().relative_to(root_str)

    return relpath


def filepath_to_objectkey(
    record: Union[File, Folder], filepath: Union[Path, UPath]
) -> str:
    root, root_str = get_storage_root_and_root_str()

    relpath = filepath_to_relpath(root=root, root_str=root_str, filepath=filepath)
    # for File, _objectkey is relative path to the storage root without suffix
    _objectkey = relpath.parent / record.name if isinstance(record, File) else relpath

    return _objectkey.as_posix()


def write_objectkey(record: sqm.SQLModel) -> None:
    """Write to _objectkey.

    An objectkey excludes the storage root and the file suffix.
    """

    def set_objectkey(record: Union[File, Folder], filepath: Union[Path, UPath]):
        _objectkey = filepath_to_objectkey(record=record, filepath=filepath)
        set_attribute(record, "_objectkey", _objectkey)

    # _local_filepath private attribute is only added
    # when creating File from data or Folder from folder
    if hasattr(record, "_local_filepath"):
        # for upsert
        if isinstance(record._local_filepath, ModelPrivateAttr):
            pass
        elif record._local_filepath is None:
            # cloud storage
            if record._cloud_filepath is not None:
                set_objectkey(record, record._cloud_filepath)
            # both _cloud_filepath and _local_filepath are None fir zarr
        # local storage
        else:
            # only set objectkey if it is configured
            if local_instance_storage_matches_local_parent(record):
                set_objectkey(record, record._local_filepath)


def upload_data_object(file) -> None:
    """Store and add file and its linked entries."""
    file_storage_key = f"{file.id}{file.suffix}"

    storage = setup_settings.instance.storage

    if file._local_filepath is not None:
        # - Look for _cloud_filepath, which is only not None if the passed filepath
        # was in the existing storage in the first place (errors within _file.py)
        # - Look for _local_filepath and check whether it's in existing storage before
        # trying to copy the file
        if (file._cloud_filepath is None) and (
            not local_instance_storage_matches_local_parent(file)
        ):
            store_object(file._local_filepath, file_storage_key)
    elif file.suffix == ".zarr" and file._memory_rep is not None:
        storagepath = storage.key_to_filepath(file_storage_key)
        print_progress = partial(print_hook, filepath=file.name)
        write_adata_zarr(file._memory_rep, storagepath, callback=print_progress)
