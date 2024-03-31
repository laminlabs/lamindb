import json
from datetime import datetime, timezone

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb._artifact import Artifact


# tested in lamin-spatial
# can't type vitessce_config because can't assume it's installed
def save_vitessce_config(vitessce_config, description: str) -> Artifact:
    """Takes a ``VitessceConfig`` object and saves it as an artifact.

    Args:
        vitessce_config (``VitessceConfig``): A VitessceConfig object.
        description: A description for the artifact.
    """
    from vitessce import VitessceConfig

    assert isinstance(vitessce_config, VitessceConfig)
    timestamp = datetime.now(timezone.utc).isoformat().split(".")[0]
    vitesse_export = f"./vitessce_export_{timestamp}.vitessce"
    vitessce_config.export(to="files", base_url="", out_dir=vitesse_export)
    logger.important(f"local export: {vitesse_export}")
    artifact = Artifact(vitesse_export, description=description)
    artifact.save()
    config_dict = vitessce_config.to_dict(base_url=artifact.path.to_url())
    logger.important(f"base url: {artifact.path.to_url()}")
    config_filename = "vitessce_config.json"
    config_file_local_path = f"{vitesse_export}/{config_filename}"
    with open(config_file_local_path, "w") as file:
        json.dump(config_dict, file)
    config_file_path = artifact.path / config_filename
    config_file_path.upload_from(config_file_local_path)
    logger.important(f"config url: {config_file_path.to_url()}")
    slug = ln_setup.settings.instance.slug
    logger.important(f"go to: https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact
