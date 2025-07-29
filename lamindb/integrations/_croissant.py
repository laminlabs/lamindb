from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

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
            raise ValueError(
                f"No valid file path found in croissant distribution: {dist}"
            )
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
