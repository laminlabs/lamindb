"""PyTorch Lightning integrations.

.. autosummary::
    :toctree: .

    LightningCallback
"""

from collections.abc import Sequence
from pathlib import Path

import lightning as pl
from lightning.pytorch import LightningModule, Trainer

import lamindb as ln


class LightningCallback(pl.Callback):
    """Saves PyTorch Lightning model checkpoints to LaminDB after each training epoch.

    Creates version families of artifacts for given `key` (relative file path).

    Args:
        path: Path to the checkpoint
        key: Artifact key in LaminDB storage
        annotate_by: Metric names to annotate artifacts with. If None, no annotation.
    """

    def __init__(
        self, path: str | Path, key: str, annotate_by: Sequence[str] | None = None
    ):
        self.path = Path(path)
        self.key = key
        self.params = annotate_by

    def on_train_start(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Validates that LaminDB Features exist for all specified params."""
        if self.params:
            missing = [
                p
                for p in self.params
                if ln.Feature.filter(name=p).one_or_none() is None
            ]
            if missing:
                s = "s" if len(missing) > 1 else ""
                raise ValueError(
                    f"Feature{s} {', '.join(missing)} missing. Create {'them' if len(missing) > 1 else 'it'} first."
                )

    def on_train_epoch_end(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Saves model checkpoint to LaminDB at end of each epoch."""
        trainer.save_checkpoint(self.path)

        af = ln.Artifact(self.path, key=self.key, kind="model").save()

        if self.params:
            values = {}
            for param in self.params:
                if param in trainer.callback_metrics:
                    v = trainer.callback_metrics[param]
                    values[param] = v.item() if hasattr(v, "item") else float(v)

            af.params.add_values(values)
            af.save()


__all__ = ["LightningCallback"]
