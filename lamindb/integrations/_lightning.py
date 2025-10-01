"""PyTorch Lightning integrations.

.. autosummary::
   :toctree: .

   LightningCallback
"""

from pathlib import Path

import lightning as pl
from lightning.pytorch import LightningModule, Trainer

import lamindb as ln


class LightningCallback(pl.Callback):
    """Saves PyTorch Lightning model checkpoints to LaminDB after each training epoch.

    Creates version families of artifacts for given `key` (relative file path).

    Args:
        key_prefix: Prefix for artifact keys in LaminDB storage
        path: Path to the checkpoint
    """

    def __init__(self, path: str | Path, key: str):
        self.path = Path(path)
        self.key = key

    def on_train_epoch_end(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Saves model checkpoint to LaminDB at end of each epoch."""
        trainer.save_checkpoint(self.path)

        ln.Artifact(
            self.path,
            key=self.key,
            kind="model",
        ).save()


__all__ = ["LightningCallback"]
