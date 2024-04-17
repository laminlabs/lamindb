from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb._artifact import Artifact

if TYPE_CHECKING:
    from vitessce import VitessceConfig


# tested & context in https://github.com/laminlabs/lamin-spatial
def save_vitessce_config(vitessce_config: VitessceConfig, description: str) -> Artifact:
    """Takes a ``VitessceConfig`` object and saves it as an artifact.

    Args:
        vitessce_config (``VitessceConfig``): A VitessceConfig object.
        description: A description for the artifact.
    """
    # can't assume vitessce is installed
    from vitessce import VitessceConfig

    # create a local _data export_ in a folder
    timestamp = datetime.now(timezone.utc).isoformat().split(".")[0]
    vitesse_export_folder = f"./vitessce_export_{timestamp}.vitessce"
    vitessce_config.export(to="files", base_url="", out_dir=vitesse_export_folder)
    logger.important(f"local export: {vitesse_export_folder}")
    # create an artifact and store the local export in th cloud
    artifact = Artifact(vitesse_export_folder, description=description)
    artifact.save()
    # create a JSON export that points to the data in the cloud
    config_dict = vitessce_config.to_dict(base_url=artifact.path.to_url())
    logger.important(f"base url: {artifact.path.to_url()}")
    # manually place that JSON export into the local data export folder
    config_filename = "vitessce_config.json"
    config_file_local_path = f"{vitesse_export_folder}/{config_filename}"
    with open(config_file_local_path, "w") as file:
        json.dump(config_dict, file)
    # manually place that JSON export into the previously registered artifact folder
    config_file_path = artifact.path / config_filename
    config_file_path.upload_from(config_file_local_path)
    # log the the URLs
    logger.important(f"config url: {config_file_path.to_url()}")
    slug = ln_setup.settings.instance.slug
    logger.important(f"go to: https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact
