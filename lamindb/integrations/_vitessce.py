from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb.models.artifact import Artifact
from lamindb.models.collection import Collection
from lamindb.models.run import Run
from lamindb.models.transform import Transform

if TYPE_CHECKING:
    from vitessce import VitessceConfig


# "unit test": https://github.com/laminlabs/lamindb/blob/main/docs/storage/vitessce.ipynb
# integration test & context: https://github.com/laminlabs/lamin-spatial/blob/main/docs/vitessce.ipynb
def save_vitessce_config(
    vitessce_config: VitessceConfig,
    key: str | None = None,
    description: str | None = None,
) -> Artifact:
    """Validates and saves a `VitessceConfig` object.

    If the `VitessceConfig` object references multiple artifacts, automatically
    creates a `Collection` and displays the "Vitessce button" next to it.

    The `VitessceConfig` artifact has `.suffix = ".vitessce.json"` and `.kind = "__lamindb_config__"`,
    which is by default hidden on the hub UI.

    Guide: :doc:`docs:vitessce`.

    Args:
        vitessce_config: A `VitessceConfig` object.
        key: A `key` for the `VitessceConfig` artifact.
        description: A `description` for the `VitessceConfig` aritifact. Is additionally
            used as `key` for a `Collection` in case the `VitessceConfig` object
            references multiple artifacts.
    """
    # can only import here because vitessce is not a dependency
    from vitessce import VitessceConfig

    assert isinstance(vitessce_config, VitessceConfig)  # noqa: S101
    vc_dict = vitessce_config.to_dict()
    try:
        url_to_artifact_dict = vitessce_config.get_artifacts()
    except AttributeError as e:
        raise SystemExit(
            "save_vitessce_config() requires vitessce>=3.4.0: pip install vitessce>=3.4.0"
        ) from e
    dataset_artifacts = list(url_to_artifact_dict.values())
    message = "\n".join([artifact.__repr__() for artifact in dataset_artifacts])
    logger.important(f"VitessceConfig references these artifacts:\n{message}")
    assert len(dataset_artifacts) > 0  # noqa: S101

    # the below will be replaced with a `ln.tracked()` decorator soon
    transform = Transform(  # type: ignore
        uid="kup03MJBsIVa0002",
        key="save_vitessce_config",
        type="function",
        version="3",
    ).save()
    run = Run(transform=transform).save()
    run.input_artifacts.set(dataset_artifacts)
    collection = None
    if len(dataset_artifacts) > 1:
        # if we have more datasets, we should create a collection
        # and attach an action to the collection
        # consicious use of description for key, see here
        # https://github.com/laminlabs/lamindb/pull/2997
        collection = Collection(dataset_artifacts, key=description).save()

    # create a JSON export
    config_file_local_path = ln_setup.settings.cache_dir / "config.vitessce.json"
    with open(config_file_local_path, "w") as file:
        json.dump(vc_dict, file)
    vitessce_config_artifact = Artifact(
        config_file_local_path,
        key=key,
        description=description,
        run=run,
        kind="__lamindb_config__",
    ).save()
    slug = ln_setup.settings.instance.slug
    logger.important(
        f"VitessceConfig: https://lamin.ai/{slug}/artifact/{vitessce_config_artifact.uid}"
    )
    if collection is None:
        # we have one and only one dataset artifact, hence the following line is OK
        dataset_artifacts[0]._actions.add(vitessce_config_artifact)
        logger.important(
            f"Dataset: https://lamin.ai/{slug}/artifact/{dataset_artifacts[0].uid}"
        )
    else:
        collection._actions.add(vitessce_config_artifact)
        logger.important(
            f"Collection: https://lamin.ai/{slug}/collection/{collection.uid}"
        )
    run.finished_at = datetime.now(timezone.utc)
    run.save()
    return vitessce_config_artifact
