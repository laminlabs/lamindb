from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb._artifact import Artifact
from lamindb._collection import Collection
from lamindb._run import Run
from lamindb._transform import Transform

if TYPE_CHECKING:
    from vitessce import VitessceConfig


# "unit test": https://github.com/laminlabs/lamindb/blob/main/docs/storage/vitessce.ipynb
# integration test & context: https://github.com/laminlabs/lamin-spatial/blob/main/docs/vitessce.ipynb
def save_vitessce_config(
    vitessce_config: VitessceConfig, description: str | None = None
) -> Artifact:
    """Validates and saves a ``VitessceConfig`` object.

    Guide: :doc:`docs:vitessce`.

    Args:
        vitessce_config (``VitessceConfig``): A `VitessceConfig` object.
        description: A description for the `VitessceConfig` artifact.

    .. versionchanged:: 0.75.1
        Now displays the "Vitessce button" on the hub next to the dataset. It additionally keeps displaying it next to the configuration file.
    .. versionchanged:: 0.70.2
        No longer saves the dataset. It only saves the `VitessceConfig` object.
    """
    # can only import here because vitessce is not a dependency
    from vitessce import VitessceConfig

    from lamindb.core.storage import VALID_SUFFIXES

    assert isinstance(vitessce_config, VitessceConfig)  # noqa: S101
    vc_dict = vitessce_config.to_dict()
    try:
        url_to_artifact_dict = vitessce_config.get_artifacts()
    except AttributeError as e:
        logger.error("Artifact registration requires vitessce package version 3.4.0 or higher.")
        raise e
    dataset_artifacts = list(url_to_artifact_dict.values())
    if len(dataset_artifacts) == 0:
        logger.warning("No artifacts were registered in this config. If intending to visualize data from artifacts, use _artifact parameters of Vitessce wrapper class constructors to facilitate registration.")
    
    # the below will be replaced with a `ln.tracked()` decorator soon
    with logger.mute():
        transform = Transform(
            uid="kup03MJBsIVa0001",
            name="save_vitessce_config",
            type="function",
            version="2",
        ).save()
    run = Run(transform=transform).save()
    run.input_artifacts.set(dataset_artifacts)
    if len(dataset_artifacts) > 1:
        # if we have more datasets, we should create a collection
        # and attach an action to the collection
        collection_of_artifacts = Collection(dataset_artifacts)
    else:
        collection_of_artifacts = None
    
    # create a JSON export
    config_file_local_path = ln_setup.settings.cache_dir / "config.vitessce.json"
    with open(config_file_local_path, "w") as file:
        json.dump(vc_dict, file)
    vitessce_config_artifact = Artifact(
        config_file_local_path, description=description, run=run
    ).save()
    if collection_of_artifacts is None:
        # we have one and only one dataset artifact, hence the following line is OK
        dataset_artifacts[0]._actions.add(vitessce_config_artifact)
    else:
        collection_of_artifacts._actions.add(vitessce_config_artifact)
    slug = ln_setup.settings.instance.slug
    logger.important(
        f"go to: https://lamin.ai/{slug}/artifact/{vitessce_config_artifact.uid}"
    )
    run.finished_at = datetime.now(timezone.utc)
    run.save()
    return vitessce_config_artifact
