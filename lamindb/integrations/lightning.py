"""PyTorch Lightning.

.. autoclass:: Callback
"""

from pathlib import Path
from typing import Any

import lightning as pl
from lightning.pytorch import LightningModule, Trainer

import lamindb as ln


class Callback(pl.Callback):
    """Saves PyTorch Lightning model checkpoints to the LaminDB instance after each training epoch.

    Creates version families of artifacts for given `key` (relative file path).

    See also: :doc:`docs:mlflow` & :doc:`docs:wandb`.

    Args:
        path: A local path to the checkpoint.
        key: The `key` for the checkpoint artifact.
        features: Features to annotate the checkpoint.

    Examples:

        Create a callback that creates artifacts for checkpoints and annotates them by the MLflow run ID::

            import lightning as pl
            from lamindb.integrations import lightning as ll

            lamindb_callback = ll.Callback(
                path=checkpoint_filename, key=artifact_key, features={"mlflow_run_id": mlflow_run.info.run_id}
            )
            trainer = pl.Trainer(callbacks=[lamindb_callback])
    """

    def __init__(
        self,
        path: str | Path,
        key: str,
        features: dict[str, Any] | None = None,
    ):
        self.path = Path(path)
        self.key = key
        self.features = features or {}

    def on_train_start(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Validates that features exist for all specified params."""
        missing = [
            feature
            for feature in self.features.keys()
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

        feature_values = dict(self.features)

        for name in self.features.keys():
            if hasattr(trainer, name):
                feature_values[name] = getattr(trainer, name)
            elif name in trainer.callback_metrics:
                metric_value = trainer.callback_metrics[name]
                feature_values[name] = (
                    metric_value.item()
                    if hasattr(metric_value, "item")
                    else float(metric_value)
                )

        if feature_values:
            af.features.add_values(feature_values)

        af.save()


__all__ = ["Callback"]
