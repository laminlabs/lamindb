from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from anndata import AnnData
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core._settings_storage import get_storage_region
from lamindb_setup.core.upath import create_path
from lnschema_core import Artifact, Run

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from tiledbsoma.io import ExperimentAmbientLabelMapping
    from upath import UPath


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


def _prepare_adatas(
    adatas: list[AnnData | UPathStr], run_uid: str | None
) -> list[AnnData]:
    add_run_uid = run_uid is not None
    adata_objects = []
    for adata in adatas:
        if isinstance(adata, AnnData):
            if add_run_uid:
                if adata.is_view:
                    raise ValueError(
                        "Can not write an `AnnData` view, please do `adata.copy()` before passing."
                    )
                else:
                    logger.warning("Mutating in-memory AnnData.")
                    adata.obs["lamin_run_uid"] = run_uid
        else:
            adata = _read_adata_h5ad_zarr(create_path(adata))
            if add_run_uid:
                adata.obs["lamin_run_uid"] = run_uid
        adata_objects.append(adata)
    return adata_objects


def save_tiledbsoma_experiment(
    adatas: list[AnnData | UPathStr],
    measurement_name: str,
    revises: Artifact | None = None,
    run: Run | None = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    append_obsm_varm: bool = False,
    artifact_kwargs: dict | None = None,
    **kwargs,
) -> Artifact:
    try:
        import tiledbsoma as soma
        import tiledbsoma.io as soma_io
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    from lamindb.core._data import get_run
    from lamindb.core.storage.paths import auto_storage_key_from_artifact_uid
    from lamindb.core.versioning import create_uid

    run = get_run(run)

    if artifact_kwargs is None:
        artifact_kwargs = {}

    appending = revises is not None

    if appending:
        _uid = None
        storepath = revises.path
    else:
        _uid, _ = create_uid(n_full_id=20)
        storage_key = auto_storage_key_from_artifact_uid(
            _uid, ".tiledbsoma", is_dir=True
        )
        storepath = setup_settings.storage.root / storage_key

    if storepath.protocol == "s3":
        ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(storepath))
    else:
        ctx = None

    storepath = storepath.as_posix()

    add_run_uid = True
    if appending:
        with soma.Experiment.open(storepath, mode="r", context=ctx) as store:
            add_run_uid = "lamin_run_uid" in store["obs"].schema.names

    adata_objects = _prepare_adatas(adatas, run.uid if add_run_uid else None)

    if appending or len(adata_objects) > 1:
        registration_mapping = soma_io.register_anndatas(
            experiment_uri=storepath if appending else None,
            adatas=adata_objects,
            measurement_name=measurement_name,
            obs_field_name=obs_id_name,
            var_field_name=var_id_name,
            append_obsm_varm=append_obsm_varm,
            context=ctx,
        )
    else:
        registration_mapping = None

    for adata_obj in adata_objects:
        soma_io.from_anndata(
            storepath,
            adata_obj,
            measurement_name,
            context=ctx,
            obs_id_name=obs_id_name,
            var_id_name=var_id_name,
            registration_mapping=registration_mapping,
            **kwargs,
        )

    return Artifact(
        storepath, run=run, revises=revises, _uid=_uid, **artifact_kwargs
    ).save()
