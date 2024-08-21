from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from anndata import AnnData
from lamin_utils import logger
from lamindb_setup.core._settings_storage import get_storage_region
from lamindb_setup.core.upath import create_path
from lnschema_core import Artifact, Run, Storage
from upath import UPath

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from tiledbsoma.io import ExperimentAmbientLabelMapping


def _read_adata_h5ad_zarr(objpath: UPath):
    from lamindb.core.storage.paths import read_adata_h5ad, read_adata_zarr

    if objpath.is_dir():
        adata = read_adata_zarr(objpath)
    else:
        adata = read_adata_h5ad(objpath)
    return adata


def _tiledb_config_s3(storepath: UPath) -> dict:
    region = get_storage_region(storepath)
    tiledb_config = {"vfs.s3.region": region}
    storage_options = storepath.storage_options
    if "key" in storage_options:
        tiledb_config["vfs.s3.aws_access_key_id"] = storage_options["key"]
    if "secret" in storage_options:
        tiledb_config["vfs.s3.aws_secret_access_key"] = storage_options["secret"]
    if "token" in storage_options:
        tiledb_config["vfs.s3.aws_session_token"] = storage_options["token"]

    return tiledb_config


def _open_tiledbsoma(
    storepath: UPath, mode: Literal["r", "w"] = "r"
) -> SOMACollection | SOMAExperiment:
    try:
        import tiledbsoma as soma
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    storepath_str = storepath.as_posix()
    if storepath.protocol == "s3":
        ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(storepath))
        # this is a strange bug
        # for some reason iterdir futher gives incorrect results
        # if cache is not invalidated
        # instead of obs and ms it gives ms and ms in the list of names
        storepath.fs.invalidate_cache()
    else:
        ctx = None

    soma_objects = [obj.name for obj in storepath.iterdir()]
    if "obs" in soma_objects and "ms" in soma_objects:
        SOMAType = soma.Experiment
    else:
        SOMAType = soma.Collection
    return SOMAType.open(storepath_str, mode=mode, context=ctx)


def register_for_tiledbsoma_store(
    storepath: UPathStr | None,
    adatas: list[AnnData | UPathStr],
    measurement_name: str,
    obs_field_name: str,
    var_field_name: str,
    append_obsm_varm: bool = False,
    run: Run | None = None,
) -> tuple[ExperimentAmbientLabelMapping, list[AnnData]]:
    try:
        import tiledbsoma as soma
        import tiledbsoma.io as soma_io
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    add_run_uid = True
    ctx = None
    if storepath is not None:
        storepath = create_path(storepath)
        if storepath.protocol == "s3":
            ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(storepath))
        if storepath.exists():
            with soma.Experiment.open(
                storepath.as_posix(), mode="r", context=ctx
            ) as store:
                add_run_uid = "lamin_run_uid" in store["obs"].schema.names
        storepath = storepath.as_posix()

    if add_run_uid:
        from lamindb.core._data import get_run

        run = get_run(run)

    adata_objects = []
    for adata in adatas:
        if isinstance(adata, AnnData):
            if add_run_uid:
                if adata.is_view:
                    raise ValueError(
                        "Can not register an `AnnData` view, please do `adata.copy()` before passing."
                    )
                else:
                    logger.warning("Mutating in-memory AnnData.")
                    adata.obs["lamin_run_uid"] = run.uid
        else:
            adata = _read_adata_h5ad_zarr(create_path(adata))
            if add_run_uid:
                adata.obs["lamin_run_uid"] = run.uid
        adata_objects.append(adata)

    registration_mapping = soma_io.register_anndatas(
        experiment_uri=storepath,
        adatas=adata_objects,
        measurement_name=measurement_name,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        append_obsm_varm=append_obsm_varm,
        context=ctx,
    )

    return registration_mapping, adata_objects


def write_tiledbsoma_store(
    store: Artifact | UPathStr,
    adata: AnnData | UPathStr,
    run: Run | None = None,
    artifact_kwargs: dict | None = None,
    **kwargs,
) -> Artifact:
    """Write `AnnData` to `tiledbsoma.Experiment`.

    Reads `AnnData`, writes it to `tiledbsoma.Experiment` and creates `lamindb.Artifact`.

    See `tiledbsoma.io.from_h5ad
    <https://tiledbsoma.readthedocs.io/en/latest/_autosummary/tiledbsoma.io.from_h5ad.html>`__.
    """
    try:
        import tiledbsoma as soma
        import tiledbsoma.io as soma_io
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    from lamindb.core._data import get_run

    if artifact_kwargs is None:
        artifact_kwargs = {}

    appending: bool = kwargs.get("registration_mapping", None) is not None
    store_is_artifact: bool = isinstance(store, Artifact)
    if store_is_artifact:
        if not appending:
            raise ValueError(
                "Trying to append to an existing store without `registration_mapping`."
            )
        storepath = store.path
    else:
        storepath = create_path(storepath)
    add_run_uid: bool = not appending

    if not isinstance(adata, AnnData):
        # create_path is used
        # in case adata is somewhere in our managed s3 bucket or just in s3
        adata = _read_adata_h5ad_zarr(create_path(adata))
    elif add_run_uid and adata.is_view:
        raise ValueError(
            "Can not write from an `AnnData` view, please do `adata.copy()` before passing."
        )

    run = get_run(run)

    if add_run_uid:
        adata.obs["lamin_run_uid"] = run.uid

    if storepath.protocol == "s3":
        ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(storepath))
    else:
        ctx = None

    soma_io.from_anndata(storepath.as_posix(), adata, context=ctx, **kwargs)

    if add_run_uid:
        del adata.obs["lamin_run_uid"]

    is_new_version_of = None
    if appending:
        if store_is_artifact:
            is_new_version_of = store
        else:
            from lamindb._artifact import (
                check_path_in_existing_storage,
                get_relative_path_to_directory,
            )

            storage = check_path_in_existing_storage(storepath)
            if isinstance(storage, Storage):
                search_by_key = get_relative_path_to_directory(
                    path=storepath, directory=UPath(storage.root)
                ).as_posix()
                is_new_version_of = Artifact.filter(
                    key=search_by_key, _key_is_virtual=False
                ).one_or_none()
                if is_new_version_of is not None:
                    logger.info(f"Assuming it is a new version of {is_new_version_of}.")

    if is_new_version_of is None:
        return Artifact(storepath, run=run, **artifact_kwargs)
    else:
        return Artifact(
            storepath, run=run, is_new_version_of=is_new_version_of, **artifact_kwargs
        )
