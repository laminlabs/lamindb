"""PyTorch Lightning integrations.

.. autosummary::
    :toctree: .

    Callback
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import lightning as pl
from lightning.pytorch import LightningModule, Trainer

import lamindb as ln


class Callback(pl.Callback):
    """Saves PyTorch Lightning model checkpoints to LaminDB after each training epoch.

    Creates version families of artifacts for given `key` (relative file path).

    Args:
        path: Path to the checkpoint
        key: Artifact key in LaminDB storage
        annotate_by: Lightning metric names to annotate artifacts with.
            Examples are 'train_loss' and 'val_loss'.
        feature_values: Additional feature values that every checkpoint gets annotated by.
            Examples are { "mlflow_run_id": mlflow_run.info.run_id }.
    """

    def __init__(
        self,
        path: str | Path,
        key: str,
        annotate_by: Sequence[str] | None = None,
        feature_values: dict[str, Any] | None = None,
    ):
        self.path = Path(path)
        self.key = key
        self.annotate_by = annotate_by
        self.feature_values = feature_values or {}

    def on_train_start(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Validates that LaminDB Features exist for all specified params."""
        all_features = list(self.annotate_by or []) + list(self.feature_values.keys())
        if all_features:
            missing = [
                feature
                for feature in all_features
                if ln.Feature.filter(name=feature).one_or_none() is None
            ]
            if missing:
                s = "s" if len(missing) > 1 else ""
                raise ValueError(
                    f"Feature{s} {', '.join(missing)} missing. Create {'them' if len(missing) > 1 else 'it'} first."
                )

    def on_train_epoch_end(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Saves model checkpoint artifacts at the end of each epoch and optionally annotates them."""
        trainer.save_checkpoint(self.path)
        af = ln.Artifact(self.path, key=self.key, kind="model").save()

        feature_values = dict(self.feature_values)

        if self.annotate_by:
            for metric_name in self.annotate_by:
                if hasattr(trainer, metric_name):
                    feature_values[metric_name] = getattr(trainer, metric_name)
                elif metric_name in trainer.callback_metrics:
                    metric_value = trainer.callback_metrics[metric_name]
                    feature_values[metric_name] = (
                        metric_value.item()
                        if hasattr(metric_value, "item")
                        else float(metric_value)
                    )

        if feature_values:
            af.features.add_values(feature_values)

        af.save()


__all__ = ["Callback"]
