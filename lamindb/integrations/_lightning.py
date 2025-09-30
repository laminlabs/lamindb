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

    Args:
        run_id: Unique identifier for this training run
        key_prefix: Prefix for artifact keys in LaminDB storage
        checkpoint_dir: Directory to save temporary checkpoint files
    """

    def __init__(
        self,
        run_id: str,
        key_prefix: str = "testmodels",
        checkpoint_dir: str | Path = "model_checkpoints",
    ):
        self.run_id = run_id
        self.key_prefix = key_prefix
        self.checkpoint_dir = Path(checkpoint_dir)

    def on_train_epoch_end(self, trainer: Trainer, pl_module: LightningModule) -> None:
        """Saves model checkpoint to LaminDB at end of each epoch."""
        checkpoint_path = (
            self.checkpoint_dir / f"{self.run_id}_epoch{trainer.current_epoch}.ckpt"
        )
        trainer.save_checkpoint(checkpoint_path)

        ln.Artifact(
            checkpoint_path,
            key=f"{self.key_prefix}/{self.run_id}.ckpt",
            kind="model",
        ).save()
