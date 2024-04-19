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

    .. versionchanged:: 0.70.2
        This function no longer saves the dataset. It only saves the VitessceConfig object.
    """
    # create a JSON export that points to the data in the cloud
    config_file_local_path = (
        ln_setup.settings.storage.cache_dir / "config.vitessce.json"
    )
    with open(config_file_local_path, "w") as file:
        json.dump(vitessce_config.to_dict(), file)
    artifact = Artifact(config_file_local_path, description=description)
    artifact.save()
    slug = ln_setup.settings.instance.slug
    logger.important(f"go to: https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact
