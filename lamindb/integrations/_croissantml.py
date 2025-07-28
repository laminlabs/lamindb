from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import lamindb as ln


def curate_from_croissantml(
    croissant_data: str | Path | dict[str, Any],
) -> ln.Artifact | ln.Collection:
    """Create annotated artifacts from a CroissantML file.

    Returns a collection if multiple files are found in `croissant_data`, otherwise a single artifact.

    Args:
        croissant_data: Path to CroissantML JSON file, JSON string, or dict

    Example:

        ::

            artifact = ln.integrations.curate_from_croissantml("dataset_metadata.json")
            collection = ln.integrations.curate_from_croissantml(croissant_dict)
    """
    import lamindb as ln

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
    description = data.get("description", "")
    version = data.get("version", "1.0")
    license_info = data.get("license", "")
    project_name = data.get("cr:projectName", "")

    # Create license feature and label if license info exists
    license_label = None
    if license_info:
        license_ulabel_type = ln.ULabel.filter(name="License", is_type=True).first()
        if not license_ulabel_type:
            license_ulabel_type = ln.ULabel(name="License", is_type=True).save()
        license_ulabel = ln.ULabel.filter(name=license_info).first()
        if not license_ulabel:
            license_label = ln.ULabel(
                name=license_info,
                description="Dataset license",
                type=license_ulabel_type,
            ).save()
    project_label = None
    if project_name:
        project_label = ln.Project.filter(name=project_name).first()
        if not project_label:
            project_label = ln.Project(name=project_name).save()

    # Extract file distributions
    artifacts = []
    file_distributions = data.get("distribution", [])
    if isinstance(file_distributions, dict):
        file_distributions = [file_distributions]
    for dist in file_distributions:
        file_id = dist.get("@id", "")
        if Path(file_id).exists():
            file_path = file_id
        else:
            content_url = dist.get("contentUrl", "")
            file_path = content_url or data.get("url", "")
        if not file_path:
            raise ValueError(f"No valid file path found in CroissantML: {dist}")
        if len(file_distributions) == 1:
            artifact_description = f"{dataset_name}"
            if file_id != dataset_name:
                artifact_description += f" ({file_id})"
            artifact_description += f" - {description}"
        else:
            artifact_description = f"{file_id}"
        artifact = ln.Artifact(  # type: ignore
            file_path,
            description=artifact_description,
            version=version,
            kind="dataset",
            run=False,
        ).save()
        if license_label:
            artifact.ulabels.add(license_label)
        if project_label:
            artifact.projects.add(project_label)
        artifacts.append(artifact)

    if len(artifacts) == 1:
        return artifacts[0]
    elif len(artifacts) > 1:
        collection = ln.Collection(  # type: ignore
            artifacts, key=dataset_name, description=description, version=version
        ).save()
        if license_label:
            collection.ulabels.add(license_label)
        if project_label:
            artifact.projects.add(project_label)
        return collection
    else:
        raise ValueError("No valid file distributions found in CroissantML data")
