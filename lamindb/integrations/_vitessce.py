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
    dataset_artifacts = []
    for vitessce_dataset in vc_dict["datasets"]:
        if "files" not in vitessce_dataset:
            raise ValueError("Each vitessce_dataset must have a 'files' key.")
        for file in vitessce_dataset["files"]:
            if "url" not in file:
                raise ValueError("Each file must have a 'url' key.")
            s3_path_last_element = file["url"].split("/")[-1]
            if not s3_path_last_element.endswith(
                (".anndata.zarr", ".zarr", ".ome.zarr")
            ):
                logger.warning(
                    "filename should end with '.anndata.zarr', '.zarr', or '.ome.zarr'."
                )
            stem_uid = (
                s3_path_last_element.replace(".anndata.zarr", "")
                .replace(".ome.zarr", "")
                .replace(".zarr", "")  # needs to come last
            )
            assert "." not in stem_uid  # noqa: S101 successfully stripped suffix
            artifact = Artifact.filter(uid__startswith=stem_uid).one_or_none()
            if artifact is None:
                logger.warning(
                    f"could not find dataset '{stem_uid}' in lamindb: {vitessce_dataset}"
                )
            else:
                dataset_artifacts.append(artifact)
    # link inputs
    with logger.mute():
        transform = Transform(name="save_vitessce_config", type="function", version="1")
        transform.save()
    run = Run(transform=transform)
    run.save()
    if len(dataset_artifacts) > 1:
        # if we have more datasets, we should create a collection
        # and attach an action to the collection
        raise NotImplementedError
    run.dataset_artifacts.set(dataset_artifacts)
    # create a JSON export
    config_file_local_path = (
        ln_setup.settings.storage.cache_dir / "config.vitessce.json"
    )
    with open(config_file_local_path, "w") as file:
        json.dump(vc_dict, file)
    vitessce_config_artifact = Artifact(
        config_file_local_path, description=description, run=run
    ).save()
    # we have one and only one dataset artifact, hence the following line is OK
    dataset_artifacts[0]._actions.add(vitessce_config_artifact)
    slug = ln_setup.settings.instance.slug
    logger.important(
        f"go to: https://lamin.ai/{slug}/artifact/{vitessce_config_artifact.uid}"
    )
    return vitessce_config_artifact
