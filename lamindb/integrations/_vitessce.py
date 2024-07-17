from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb._artifact import Artifact
from lamindb._run import Run
from lamindb._transform import Transform

if TYPE_CHECKING:
    from vitessce import VitessceConfig


# tested & context in https://github.com/laminlabs/lamin-spatial
def save_vitessce_config(vitessce_config: VitessceConfig, description: str) -> Artifact:
    """Validates and saves a ``VitessceConfig`` object.

    Example: :doc:`docs:vitessce`.

    Args:
        vitessce_config (``VitessceConfig``): A VitessceConfig object.
        description: A description for the artifact.

    .. versionchanged:: 0.70.2
        This function no longer saves the dataset. It only saves the VitessceConfig object.
    """
    vc_dict = vitessce_config.to_dict()
    # validate
    datasets = vc_dict["datasets"]
    input_artifacts = []
    for dataset in datasets:
        if "files" not in dataset:
            raise ValueError("Each dataset must have a 'files' key.")
        for file in dataset["files"]:
            if "url" not in file:
                raise ValueError("Each file must have a 'url' key.")
            filename = file["url"].split("/")[-1]
            if not filename.endswith(
                (".anndata.zarr", ".spatialdata.zarr", ".ome.zarr")
            ):
                raise ValueError(
                    "filename must end with '.anndata.zarr', '.spatialdata.zarr', or '.ome.zarr'."
                )
            filestem = (
                filename.replace(".anndata.zarr", "")
                .replace(".spatialdata.zarr", "")
                .replace(".ome.zarr", "")
            )
            artifact = Artifact.filter(uid__startswith=filestem).one_or_none()
            if artifact is None:
                logger.warning(
                    f"could not find dataset '{filestem}' in lamindb: {dataset}"
                )
            else:
                input_artifacts.append(artifact)
    # link inputs
    with logger.mute():
        transform = Transform(name="save_vitessce_config", type="function", version="1")
        transform.save()
    run = Run(transform=transform)
    run.save()
    run.input_artifacts.set(input_artifacts)
    # create a JSON export
    config_file_local_path = (
        ln_setup.settings.storage.cache_dir / "config.vitessce.json"
    )
    with open(config_file_local_path, "w") as file:
        json.dump(vc_dict, file)
    artifact = Artifact(config_file_local_path, description=description, run=run)
    artifact.save()
    slug = ln_setup.settings.instance.slug
    logger.important(f"go to: https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact
