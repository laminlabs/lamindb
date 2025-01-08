from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pandas as pd
import pyarrow as pa
from anndata import AnnData, read_h5ad
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core._settings_storage import get_storage_region
from lamindb_setup.core.upath import LocalPathClasses, create_path
from packaging import version

from lamindb.models import Artifact, Run

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from upath import UPath


def _load_h5ad_zarr(objpath: UPath):
    from lamindb.core.loaders import load_anndata_zarr, load_h5ad

    if objpath.is_dir():
        adata = load_anndata_zarr(objpath)
    else:
        # read only local in backed for now
        # in principle possible to read remote in backed also
        if isinstance(objpath, LocalPathClasses):
            adata = read_h5ad(objpath.as_posix(), backed="r")
        else:
            adata = load_h5ad(objpath)
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


def save_tiledbsoma_experiment(
    # Artifact args
    adatas: list[AnnData | UPathStr],
    key: str | None = None,
    description: str | None = None,
    run: Run | None = None,
    revises: Artifact | None = None,
    # tiledbsoma.io.from_anndata args
    measurement_name: str = "RNA",
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    append_obsm_varm: bool = False,
    # additional keyword args for tiledbsoma.io.from_anndata
    **kwargs,
) -> Artifact:
    """Write `AnnData` to `tiledbsoma.Experiment`.

    Reads `AnnData` objects, writes them to `tiledbsoma.Experiment`, creates & saves an {class}`~lamindb.Artifact`.

    Populates a column `lamin_run_uid` column in `obs` with the current `run.uid`.

    Is based on `tiledbsoma.io.from_anndata
    <https://tiledbsoma.readthedocs.io/en/latest/_autosummary/tiledbsoma.io.from_anndata.html>`__.

    Args:
        adatas: `AnnData` objects to write, in-memory or on-disk.
        key: An optional key to reference the artifact.
        description: A description.
        run: The run that creates the artifact.
        revises: `lamindb.Artifact` with `tiledbsoma.Experiment` to append to.
        measurement_name: The name of the measurement to store data in `tiledbsoma.Experiment`.
        obs_id_name: Which `AnnData` `obs` column to use for append mode.
        var_id_name: Which `AnnData` `var` column to use for append mode.
        append_obsm_varm: Whether to append `obsm` and `varm` in append mode .
        **kwargs: Keyword arguments passed to `tiledbsoma.io.from_anndata`.
    """
    try:
        import tiledbsoma as soma
        import tiledbsoma.io as soma_io
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e

    from lamindb.core._data import get_run
    from lamindb.core.storage.paths import auto_storage_key_from_artifact_uid
    from lamindb.core.versioning import create_uid

    run = get_run(run)

    appending = revises is not None
    if appending:
        storepath = revises.path
    else:
        uid, _ = create_uid(n_full_id=20)
        storage_key = auto_storage_key_from_artifact_uid(
            uid, ".tiledbsoma", is_dir=True
        )
        storepath = setup_settings.storage.root / storage_key

    if storepath.protocol == "s3":
        ctx = soma.SOMATileDBContext(tiledb_config=_tiledb_config_s3(storepath))
    else:
        ctx = None

    storepath = storepath.as_posix()

    add_run_uid = True
    run_uid_dtype = "category"
    if appending:
        with soma.Experiment.open(storepath, mode="r", context=ctx) as store:
            obs_schema = store["obs"].schema
            add_run_uid = "lamin_run_uid" in obs_schema.names
            # this is needed to enable backwards compatibility with tiledbsoma stores
            # created before PR 2300
            if add_run_uid:
                column_type = obs_schema.types[obs_schema.names.index("lamin_run_uid")]
                if not isinstance(column_type, pa.DictionaryType):
                    run_uid_dtype = None

    if add_run_uid and run is None:
        raise ValueError("Pass `run`")

    adata_objects = []
    for adata in adatas:
        if isinstance(adata, AnnData):
            if add_run_uid and adata.is_view:
                raise ValueError(
                    "Can not write an `AnnData` view, please do `adata.copy()` before passing."
                )
        else:
            adata = _load_h5ad_zarr(create_path(adata))
        if add_run_uid:
            adata.obs["lamin_run_uid"] = pd.Series(
                run.uid, index=adata.obs.index, dtype=run_uid_dtype
            )
        adata_objects.append(adata)

    registration_mapping = kwargs.get("registration_mapping", None)
    if registration_mapping is None and (appending or len(adata_objects) > 1):
        registration_mapping = soma_io.register_anndatas(
            experiment_uri=storepath if appending else None,
            adatas=adata_objects,
            measurement_name=measurement_name,
            obs_field_name=obs_id_name,
            var_field_name=var_id_name,
            append_obsm_varm=append_obsm_varm,
            context=ctx,
        )

    resize_experiment = False
    if registration_mapping is not None:
        if version.parse(soma.__version__) < version.parse("1.15.0rc4"):
            n_observations = len(registration_mapping.obs_axis.data)
        else:
            n_observations = registration_mapping.get_obs_shape()
            resize_experiment = True
    else:  # happens only if not appending and only one adata passed
        assert len(adata_objects) == 1  # noqa: S101
        n_observations = adata_objects[0].n_obs

    logger.important(f"Writing the tiledbsoma store to {storepath}")
    for adata_obj in adata_objects:
        if resize_experiment and soma.Experiment.exists(storepath, context=ctx):
            # can only happen if registration_mapping is not None
            soma_io.resize_experiment(
                storepath,
                nobs=n_observations,
                nvars=registration_mapping.get_var_shapes(),
                context=ctx,
            )
            resize_experiment = False
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

    artifact = Artifact(
        storepath,
        key=key,
        description=description,
        run=run,
        revises=revises,
        _is_internal_call=True,
    )
    artifact.n_observations = n_observations
    artifact.otype = "tiledbsoma"

    return artifact.save()
