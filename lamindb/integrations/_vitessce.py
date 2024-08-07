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
    valid_composite_zarr_suffixes = [
        suffix for suffix in VALID_SUFFIXES.COMPOSITE if suffix.endswith(".zarr")
    ]
    # validate
    dataset_artifacts = []
    for vitessce_dataset in vc_dict["datasets"]:
        # didn't find good ways to violate the below, hence using plain asserts
        # without user feedback
        assert "files" in vitessce_dataset  # noqa: S101
        assert vitessce_dataset["files"]  # noqa: S101
        for file in vitessce_dataset["files"]:
            if "url" not in file:
                raise ValueError("Each file must have a 'url' key.")
            s3_path = file["url"]
            s3_path_last_element = s3_path.split("/")[-1]
            # note 1: the following parses the stem uid of the artifact from the S3 path
            # there might be a better way of doing this in case the vitessce config
            # gets updated in the future; but given these paths are set in stone
            # this should be more robust than it looks
            #
            # note 2: what's not great is the fact that people might use composite suffixes we don't recognize
            # I don't know what to do about it other than documenting it clearly
            # https://github.com/laminlabs/lamindb/blob/main/lamindb/core/storage/_valid_suffixes.py
            # https://docs.lamin.ai/lamindb.core.storage.valid_suffixes
            #
            # now start with attempting to strip the composite suffix candidates
            for suffix in valid_composite_zarr_suffixes:
                s3_path_last_element = s3_path_last_element.replace(suffix, "")
            # in case there was no hit, strip plain ".zarr"
            artifact_stem_uid = s3_path_last_element.replace(".zarr", "")
            # if there is still a "." in string, we
            if "." in artifact_stem_uid:
                raise ValueError(
                    f"Suffix should be '.zarr' or one of {valid_composite_zarr_suffixes}. Inspect your path {s3_path}"
                )
            artifact = Artifact.filter(uid__startswith=artifact_stem_uid).one_or_none()
            if artifact is None:
                raise ValueError(
                    f"Could not find dataset with stem uid '{artifact_stem_uid}' in lamindb: {vitessce_dataset}. Did you follow https://docs.lamin.ai/vitessce? It appears the AWS S3 path doesn't encode a lamindb uid."
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
    run.input_artifacts.set(dataset_artifacts)
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
