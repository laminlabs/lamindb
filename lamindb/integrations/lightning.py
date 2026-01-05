"""PyTorch Lightning integration for LaminDB.

.. autoclass:: LaminCheckpoint
.. autoclass:: SaveConfigCallback
.. autofunction:: register_features
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import torch
from lightning.fabric.utilities.cloud_io import get_filesystem
from lightning.pytorch.callbacks.model_checkpoint import ModelCheckpoint
from lightning.pytorch.cli import SaveConfigCallback as _SaveConfigCallback

import lamindb as ln

if TYPE_CHECKING:
    from datetime import timedelta

    import lightning.pytorch as pl
    from lightning.fabric.utilities.types import _PATH


def save_lightning_features() -> None:
    """Register LaminDB features used by LaminCheckpoint.

    Creates the following features if they don't already exist:
        - is_best_model (bool): Whether this checkpoint is the best model.
        - score (float): The monitored metric score.
        - model_rank (int): Rank among all checkpoints (0 = best).
        - logger_name (str): Name from the first Lightning logger.
        - logger_version (str): Version from the first Lightning logger.

    Example::

        from lamindb.integrations import lightning as ll

        ll.save_lightning_features()
    """
    ln.Feature(name="is_best_model", dtype=bool).save()
    ln.Feature(name="score", dtype=float).save()
    ln.Feature(name="model_rank", dtype=int).save()
    ln.Feature(name="logger_name", dtype=str).save()
    ln.Feature(name="logger_version", dtype=str).save()


class LaminCheckpoint(ModelCheckpoint):
    """ModelCheckpoint that uploads checkpoints to LaminDB.

    Extends Lightning's ModelCheckpoint with LaminDB artifact tracking.
    Checkpoints are versioned under a single key.

    Args:
        key: The artifact key for checkpoints in LaminDB.
        features: Features to annotate checkpoints with.
            Values can be static or None (auto-populated from trainer metrics/attributes).
        dirpath: Local directory for checkpoints.
        monitor: Quantity to monitor for saving best checkpoint.
        verbose: Verbosity mode.
        save_last: Save a copy of the last checkpoint.
        save_top_k: Number of best checkpoints to keep.
        save_weights_only: Save only model weights (not optimizer state).
        mode: One of "min" or "max" for monitor comparison.
        auto_insert_metric_name: Include metric name in checkpoint filename.
        every_n_train_steps: Checkpoint every N training steps.
        train_time_interval: Checkpoint at time intervals.
        every_n_epochs: Checkpoint every N epochs.
        save_on_train_epoch_end: Run checkpointing at end of training epoch.
        enable_version_counter: Append version to filename to avoid collisions.

    Example::

        import lightning as pl
        from lamindb.integrations import lightning as ll

        # Pre-create features
        ll.save_lightning_features()

        callback = ll.LaminCheckpoint(
            key="experiments/my_model.ckpt",
            features={
                "val_loss": None,  # auto-populated from trainer
            },
            monitor="val_loss",
            save_top_k=3,
            mode="min",
        )
        trainer = pl.Trainer(callbacks=[callback])
    """

    def __init__(
        self,
        key: str,
        *,
        features: dict[str, Any] | None = None,
        dirpath: _PATH | None = None,
        monitor: str | None = None,
        verbose: bool = False,
        save_last: bool | None = None,
        save_top_k: int = 1,
        save_weights_only: bool = False,
        mode: Literal["min", "max"] = "min",
        auto_insert_metric_name: bool = True,
        every_n_train_steps: int | None = None,
        train_time_interval: timedelta | None = None,
        every_n_epochs: int | None = None,
        save_on_train_epoch_end: bool | None = None,
        enable_version_counter: bool = True,
    ) -> None:
        super().__init__(
            dirpath=dirpath,
            monitor=monitor,
            verbose=verbose,
            save_last=save_last,
            save_top_k=save_top_k,
            save_weights_only=save_weights_only,
            mode=mode,
            auto_insert_metric_name=auto_insert_metric_name,
            every_n_train_steps=every_n_train_steps,
            train_time_interval=train_time_interval,
            every_n_epochs=every_n_epochs,
            save_on_train_epoch_end=save_on_train_epoch_end,
            enable_version_counter=enable_version_counter,
        )
        self.key = key
        self.features = features or {}

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Validate that all specified features exist."""
        super().setup(trainer, pl_module, stage)
        if trainer.is_global_zero:
            self._validate_features()

    def _validate_features(self) -> None:
        """Validate that user-specified features exist."""
        missing = [
            name
            for name in self.features
            if ln.Feature.filter(name=name).one_or_none() is None
        ]
        if missing:
            s = "s" if len(missing) > 1 else ""
            raise ValueError(
                f"Feature{s} {', '.join(missing)} missing. "
                f"Create {'them' if len(missing) > 1 else 'it'} first."
            )

    def _save_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Save checkpoint locally and upload to LaminDB."""
        super()._save_checkpoint(trainer, filepath)
        if trainer.is_global_zero:
            self._upload(trainer, filepath)

    def _upload(self, trainer: pl.Trainer, filepath: str) -> None:
        """Upload checkpoint artifact with associated features."""
        artifact = ln.Artifact(
            filepath, key=self.key, kind="model", description="Model checkpoint"
        )
        artifact.save()

        feature_values: dict[str, Any] = {}

        # Collect auto-populated values for requested features
        if "logger_name" in self.features and trainer.loggers:
            feature_values["logger_name"] = trainer.loggers[0].name
        if "logger_version" in self.features and trainer.loggers:
            version = trainer.loggers[0].version
            feature_values["logger_version"] = (
                version if isinstance(version, str) else f"version_{version}"
            )

        is_best = self.best_model_path == filepath
        if "is_best_model" in self.features:
            if is_best:
                self._clear_best_model_flags()
            feature_values["is_best_model"] = is_best

        if "score" in self.features and self.current_score is not None:
            score = self.current_score
            if torch.is_tensor(score):
                score = score.item()
            feature_values["score"] = float(score)

        # User-specified features (non-auto)
        for name, value in self.features.items():
            if name in {
                "is_best_model",
                "score",
                "model_rank",
                "logger_name",
                "logger_version",
            }:
                continue
            if value is not None:
                feature_values[name] = value
            elif hasattr(trainer, name):
                feature_values[name] = getattr(trainer, name)
            elif name in trainer.callback_metrics:
                metric = trainer.callback_metrics[name]
                feature_values[name] = (
                    metric.item() if hasattr(metric, "item") else float(metric)
                )

        if feature_values:
            artifact.features.add_values(feature_values)

        if "model_rank" in self.features:
            self._update_model_ranks()

    def _clear_best_model_flags(self) -> None:
        """Set is_best_model=False on previous best checkpoints."""
        for artifact in ln.Artifact.filter(key=self.key):
            vals = artifact.features.get_values()
            if vals.get("is_best_model"):
                artifact.features.add_values({"is_best_model": False})

    def _update_model_ranks(self) -> None:
        """Update model_rank feature for all checkpoints under this key."""
        artifacts = list(ln.Artifact.filter(key=self.key))
        scored = []
        for af in artifacts:
            vals = af.features.get_values()
            if "score" in vals:
                scored.append((vals["score"], af))

        scored.sort(key=lambda x: x[0], reverse=(self.mode == "max"))

        for rank, (_, af) in enumerate(scored):
            af.features.add_values({"model_rank": rank})


class SaveConfigCallback(_SaveConfigCallback):
    """SaveConfigCallback that also uploads config to the instance.

    Use with LightningCLI to save the resolved config.yaml alongside checkpoints.

    Example::

        from lightning.pytorch.cli import LightningCLI
        from lamindb.integrations.lightning import SaveConfigCallback

        cli = LightningCLI(
            MyModel,
            MyDataModule,
            save_config_callback=SaveConfigCallback,
        )
    """

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Save resolved config.yaml locally and upload alongside checkpoints."""
        if self.already_saved:  # type: ignore
            return

        log_dir = trainer.log_dir
        if self.save_to_log_dir and log_dir is not None:
            config_path = Path(log_dir) / self.config_filename
            fs = get_filesystem(log_dir)

            if not self.overwrite:
                file_exists = (
                    fs.isfile(config_path) if trainer.is_global_zero else False
                )
                file_exists = trainer.strategy.broadcast(file_exists)
                if file_exists:
                    raise RuntimeError(f"Config file already exists: {config_path}")

            if trainer.is_global_zero:
                fs.makedirs(log_dir, exist_ok=True)
                self.parser.save(
                    self.config,
                    config_path,
                    skip_none=False,
                    overwrite=self.overwrite,
                    multifile=self.multifile,
                )
                self._upload_config(trainer, config_path)

        if trainer.is_global_zero:
            self.save_config(trainer, pl_module, stage)
            self.already_saved = True

        self.already_saved = trainer.strategy.broadcast(self.already_saved)

    def _upload_config(self, trainer: pl.Trainer, config_path: str) -> None:
        """Upload config under same key prefix as checkpoints."""
        checkpoint_cb = self._get_checkpoint_callback(trainer)
        if checkpoint_cb is None:
            return

        key = str(Path(checkpoint_cb.key).parent / self.config_filename)
        artifact = ln.Artifact(config_path, key=key, description="Lightning CLI config")
        artifact.save()

    def _get_checkpoint_callback(self, trainer: pl.Trainer) -> LaminCheckpoint | None:
        """Find LaminCheckpoint callback if present."""
        for cb in trainer.callbacks:
            if isinstance(cb, LaminCheckpoint):
                return cb
        return None


__all__ = ["LaminCheckpoint", "SaveConfigCallback", "save_lightning_features"]
