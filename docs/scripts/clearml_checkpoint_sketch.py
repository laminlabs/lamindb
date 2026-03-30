"""Sketch for a future external `lamin-clearml` package.

This file is intentionally not imported by `lamindb`. It captures a concrete
design for an external package that layers ClearML model registration on top of
the artifact lifecycle API exposed by `lamindb.integrations.lightning`.

Suggested future package layout:

- `lamin_clearml/lightning.py`
  - `ClearMLArtifactObserver`
  - `ClearMLCheckpoint`
- `lamin_clearml/cli.py`
  - project-specific CLI glue if needed

The central design choice is composition first: an observer consumes
`ArtifactSavedEvent` and `ArtifactRemovedEvent` from Lamin. A convenience
`ClearMLCheckpoint` can then subclass Lamin's `Checkpoint` only to attach that
observer automatically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from lamindb.integrations.lightning import (
    ArtifactRemovedEvent,
    ArtifactSavedEvent,
    Checkpoint,
)

if TYPE_CHECKING:
    from clearml import OutputModel, Task


@dataclass
class ClearMLArtifactObserver:
    """Bridge Lamin artifact events to ClearML model registration.

    This observer uses the stable `storage_uri` exposed by Lamin's artifact
    events. For cloud-backed artifacts, the URI is typically an `s3://...`
    path suitable for `OutputModel.update_weights(register_uri=...)`.
    """

    task: Task
    framework: str = "PyTorch"
    model_name_prefix: str | None = None
    delete_evicted_models: bool = False
    output_models_by_key: dict[str, OutputModel] = field(default_factory=dict)

    def on_artifact_saved(self, event: ArtifactSavedEvent) -> None:
        """Register checkpoint artifacts and optionally upload config metadata."""
        if event.kind == "checkpoint":
            output_model = self._create_output_model(event)
            output_model.update_weights(
                weights_filename=event.local_path.name,
                register_uri=event.storage_uri,
                iteration=event.trainer.global_step,
                auto_delete_file=False,
                async_enable=True,
            )
            self.output_models_by_key[event.key] = output_model
            return

        if event.kind == "config":
            self.task.upload_artifact(
                name="lightning-config",
                artifact_object=event.storage_uri,
            )
            return

        if event.kind == "hparams":
            self.task.upload_artifact(
                name="lightning-hparams",
                artifact_object=event.storage_uri,
            )

    def on_artifact_removed(self, event: ArtifactRemovedEvent) -> None:
        """Optionally mirror top-k eviction in ClearML.

        Lamin removal events are currently emitted when Lightning removes the
        local checkpoint file. A future external package can decide whether that
        should also remove or archive the matching ClearML `OutputModel`.
        """
        if not self.delete_evicted_models:
            return
        if event.kind != "checkpoint":
            return

        output_model = self.output_models_by_key.pop(event.key, None)
        if output_model is None:
            return

        # Example policy hooks for a future package:
        # output_model.archive()
        # output_model.delete()

    def _create_output_model(self, event: ArtifactSavedEvent) -> OutputModel:
        from clearml import OutputModel

        return OutputModel(
            task=self.task,
            framework=self.framework,
            name=self._model_name_for_event(event),
        )

    def _model_name_for_event(self, event: ArtifactSavedEvent) -> str:
        stem = Path(event.key).stem
        if self.model_name_prefix is None:
            return stem
        return f"{self.model_name_prefix}-{stem}"


class ClearMLCheckpoint(Checkpoint):
    """Convenience subclass for projects that prefer inheritance.

    Composition remains the recommended implementation strategy, but a concrete
    subclass is useful for YAML-driven trainer configuration and for projects
    that want a single callback entry point.
    """

    def __init__(
        self,
        *,
        task: Task,
        model_name_prefix: str | None = None,
        delete_evicted_models: bool = False,
        artifact_observers: list[Any] | None = None,
        **kwargs: Any,
    ) -> None:
        observers = list(artifact_observers or [])
        observers.append(
            ClearMLArtifactObserver(
                task=task,
                model_name_prefix=model_name_prefix,
                delete_evicted_models=delete_evicted_models,
            )
        )
        super().__init__(artifact_observers=observers, **kwargs)


def example_usage(task: Task) -> dict[str, Any]:
    """Return a minimal wiring example for documentation or tests."""
    checkpoint = ClearMLCheckpoint(
        task=task,
        monitor="val_loss",
        mode="min",
        save_top_k=3,
    )

    # A future external package can continue to rely on Lamin's
    # SaveConfigCallback because config artifacts are already routed through the
    # checkpoint artifact lifecycle pipeline.
    return {
        "callbacks": [checkpoint],
        "notes": [
            "Use lamindb.integrations.lightning.SaveConfigCallback alongside this checkpoint.",
            "Checkpoint, config, and hparams artifacts will all be visible to ClearML through the observer.",
        ],
    }