from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup.core.upath import UPath

if TYPE_CHECKING:
    import lamindb as ln


def curate_from_croissant(
    croissant_data: str | Path | dict[str, Any],
    run: ln.Run | None = None,
) -> ln.Artifact | ln.Collection:
    """Create annotated artifacts from a CroissantML file.

    Returns a collection if multiple files are found in `croissant_data`, otherwise a single artifact.

    Args:
        croissant_data: Path to CroissantML JSON file or dictionary.

    Example:

        ::

            artifact = ln.integrations.curate_from_croissant("dataset_metadata.json")
    """
    import lamindb as ln

    from ..models.artifact import check_path_in_existing_storage

    # Load CroissantML data
    if isinstance(croissant_data, (str, Path)):
        if not Path(croissant_data).exists():
            raise FileNotFoundError(f"File not found: {croissant_data}")
        with open(croissant_data, encoding="utf-8") as f:
            data = json.load(f)
    elif isinstance(croissant_data, dict):
        data = croissant_data
    else:
        raise ValueError(
            "croissant_data must be a file path, JSON string, or dictionary"
        )

    # Validate basic structure
    if data.get("@type") != "Dataset":
        raise ValueError("CroissantML @type must be 'Dataset'")

    if "name" not in data:
        raise ValueError("CroissantML must have a 'name' field")

    # Extract basic metadata
    dataset_name = data["name"]
    description = data.get("description", None)
    version = data.get("version", None)
    license_info = data.get("license", None)
    project_name = data.get("cr:projectName", None)

    # Create license feature and label if license info exists
    license_label = None
    if license_info:
        license_label_type = ln.ULabel.filter(name="License", is_type=True).first()
        if not license_label_type:
            license_label_type = ln.ULabel(name="License", is_type=True).save()
        license_label = ln.ULabel.filter(name=license_info).first()
        if not license_label:
            license_label = ln.ULabel(
                name=license_info,
                description="Dataset license",
                type=license_label_type,
            ).save()
    project_label = None
    if project_name:
        project_label = ln.Project.filter(name=project_name).first()
        if not project_label:
            project_label = ln.Project(name=project_name).save()

    # Extract file distributions
    artifacts = []
    file_distributions = data.get("distribution", [])
    if not file_distributions:
        raise ValueError("No file distributions found in croissant data")
    for dist in file_distributions:
        file_id = dist.get("@id", "")
        if Path(file_id).exists():
            file_path = file_id
        else:
            content_url = dist.get("contentUrl", "")
            file_path = content_url or data.get("url", "")
        if not file_path:
            raise ValueError(f"No file path found in croissant distribution: {dist}")
        if not UPath(file_path).exists():
            raise ValueError(f"Inferred file path does not exist: {file_path}")
        result = check_path_in_existing_storage(
            file_path, check_hub_register_storage=ln_setup.settings.instance.is_on_hub
        )
        if isinstance(result, ln.Storage):
            key = None  # will automatically use existing storage key
        else:
            current_storage_location = (
                ln.settings.storage
                if not ln.setup.settings.instance.keep_artifacts_local
                else ln.settings.local_storage
            )
            logger.warning(
                f"file path {file_path} is not part of a known storage location, will be duplicated to: {current_storage_location}"
            )
            key = file_id
        if len(file_distributions) == 1:
            # it doesn't make sense to have the dataset name on the individual
            # artifact if it's part of a collection
            artifact_description = dataset_name
            if description is not None:
                artifact_description += f" - {description}"
        else:
            artifact_description = None
        artifact = ln.Artifact(  # type: ignore
            file_path,
            key=key,
            description=artifact_description,
            version=version,
            kind="dataset",
            run=run,
        ).save()
        if license_label:
            artifact.ulabels.add(license_label)
        if project_label:
            artifact.projects.add(project_label)
        artifacts.append(artifact)

    if len(artifacts) == 1:
        return artifacts[0]
    else:
        collection = ln.Collection(  # type: ignore
            artifacts, key=dataset_name, description=description, version=version
        ).save()
        if license_label:
            collection.ulabels.add(license_label)
        if project_label:
            collection.projects.add(project_label)
        return collection
