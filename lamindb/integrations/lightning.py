"""PyTorch Lightning integration for LaminDB.

.. autoclass:: Checkpoint
.. autoclass:: SaveConfigCallback
.. autofunction:: save_lightning_features
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Literal

import lightning.pytorch as pl
from lamin_utils import logger
from lightning.fabric.utilities.cloud_io import get_filesystem
from lightning.pytorch.callbacks.model_checkpoint import ModelCheckpoint
from lightning.pytorch.cli import SaveConfigCallback as _SaveConfigCallback

import lamindb as ln

if TYPE_CHECKING:
    from datetime import timedelta

    from lightning.fabric.utilities.types import _PATH


_RUN_AUTO_FEATURES: Final = frozenset({"logger_name", "logger_version"})
_ARTIFACT_AUTO_FEATURES: Final = frozenset({"is_best_model", "score", "model_rank"})
_SUPPORTED_AUTO_FEATURES: Final = _RUN_AUTO_FEATURES | _ARTIFACT_AUTO_FEATURES


def save_lightning_features() -> None:
    """Register LaminDB features used by the Lightning integration Checkpoint.

    Creates the following features if they do not already exist:

    - lamindb.lightning (feature type): Parent feature type for the below lightning features.
    - is_best_model (bool): Whether this checkpoint is the best model.
    - score (float): The monitored metric score.
    - model_rank (int): Rank among all checkpoints (0 = best).
    - logger_name (str): Name from the first Lightning logger.
    - logger_version (str): Version from the first Lightning logger.

    Example::

        from lamindb.integrations import lightning as ll

        ll.save_lightning_features()
    """
    # lamindb.lightning is lowercase and would trigger a naming warning -> mute
    with logger.mute():
        # normal matching fails because of non-matching dtype (__lamindb_lightning__ vs None)
        if (
            lightning_feature_type := ln.Feature.filter(
                name="lamindb.lightning"
            ).one_or_none()
        ) is None:
            lightning_feature_type = ln.Feature(  # type: ignore[call-overload]
                name="lamindb.lightning",
                description="Auto-generated features tracking lightning parameters & metrics",
                is_type=True,
            )
            lightning_feature_type._dtype_str = "__lamindb_lightning__"
            lightning_feature_type.save()

    ln.Feature(name="is_best_model", dtype=bool, type=lightning_feature_type).save()
    ln.Feature(name="score", dtype=float, type=lightning_feature_type).save()
    ln.Feature(name="model_rank", dtype=int, type=lightning_feature_type).save()
    ln.Feature(name="logger_name", dtype=str, type=lightning_feature_type).save()
    ln.Feature(name="logger_version", dtype=str, type=lightning_feature_type).save()


class Checkpoint(ModelCheckpoint):
    """ModelCheckpoint that annotates torch lightning checkpoints.

    Extends Lightning's ModelCheckpoint with artifact creation & feature annotation.
    Checkpoints are stored at semantic paths like `{dirpath}/epoch=0-val_loss=0.5.ckpt`.
    Each checkpoint is a separate artifact.
    Query with `ln.Artifact.filter(key__startswith=callback.dirpath)`.

    If available in the instance, the following features are automatically tracked:
    `is_best_model`, `score`, `model_rank`, `logger_name`, `logger_version`.

    Args:
        dirpath: Directory for checkpoints (reflected in cloud paths).
        features: Features to annotate checkpoints with.
            Values can be static or None (auto-populated from trainer metrics/attributes).
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

    Examples:

        Using the API::

            import lightning as pl
            from lamindb.integrations import lightning as ll

            # Optional one-time setup to enable automated lightning specific feature tracking
            ll.save_lightning_features()

            callback = ll.Checkpoint(
                dirpath="deployments/my_model/",
                monitor="val_loss",
                save_top_k=3,
            )

            trainer = pl.Trainer(callbacks=[callback])
            trainer.fit(model, dataloader)

            # Query checkpoints
            ln.Artifact.filter(key__startswith=callback.dirpath)



        Using the CLI::

            # config.yaml
            trainer:
              callbacks:
                - class_path: lamindb.integrations.lightning.Checkpoint
                  init_args:
                    dirpath: deployments/my_model/
                    monitor: val_loss
                    save_top_k: 3

            # Run with:
            # python main.py fit --config config.yaml
    """

    def __init__(
        self,
        dirpath: _PATH,
        *,
        features: dict[str, Any] | None = None,
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
        self.features = features or {}
        self._available_auto_features: set[str] = set()
        self._run_features_added = False

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Validate user features and detect available auto-features."""
        super().setup(trainer, pl_module, stage)

        if trainer.is_global_zero:
            # Validate user-specified features exist
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

            # Detect which auto-features are available
            for name in _SUPPORTED_AUTO_FEATURES:
                if ln.Feature.filter(name=name).one_or_none() is not None:
                    self._available_auto_features.add(name)

    def _get_artifact_key(self, filepath: str) -> str:
        """Return the artifact key for this checkpoint."""
        prefix = self.dirpath.rstrip("/")
        return f"{prefix}/{Path(filepath).name}"

    def _get_key_filter(self) -> dict[str, str]:
        """Return filter kwargs for querying artifacts from this callback."""
        return {"key__startswith": self.dirpath.rstrip("/") + "/"}

    def _save_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Save checkpoint to the instance."""
        super()._save_checkpoint(trainer, filepath)

        if trainer.is_global_zero:
            # Run-level features (once per run)
            if ln.context.run and not self._run_features_added:
                run_features = {}
                if "logger_name" in self._available_auto_features and trainer.loggers:
                    run_features["logger_name"] = trainer.loggers[0].name
                if (
                    "logger_version" in self._available_auto_features
                    and trainer.loggers
                ):
                    version = trainer.loggers[0].version
                    run_features["logger_version"] = (
                        version if isinstance(version, str) else f"version_{version}"
                    )
                if run_features:
                    ln.context.run.features.add_values(run_features)
                self._run_features_added = True

            artifact = ln.Artifact(
                filepath,
                key=self._get_artifact_key(filepath),
                kind="model",
                description="Lightning model checkpoint",
            )
            artifact.save()

            # Artifact-level features
            feature_values: dict[str, Any] = {}

            is_best = self.best_model_path == filepath
            if "is_best_model" in self._available_auto_features:
                if is_best:
                    self._clear_best_model_flags()
                feature_values["is_best_model"] = is_best

            # lazy import for faster import of the class
            import torch

            if (
                "score" in self._available_auto_features
                and self.current_score is not None
            ):
                score = self.current_score
                if torch.is_tensor(score):
                    score = score.item()
                feature_values["score"] = float(score)

            # User-specified features -> for now also added to artifacts
            for name, value in self.features.items():
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

            if "model_rank" in self._available_auto_features:
                self._update_model_ranks()

    def _clear_best_model_flags(self) -> None:
        """Set is_best_model=False on previous best checkpoints."""
        for artifact in ln.Artifact.filter(**self._get_key_filter()):
            vals = artifact.features.get_values()
            if vals.get("is_best_model"):
                artifact.features.remove_values("is_best_model", value=True)
                artifact.features.add_values({"is_best_model": False})

    def _update_model_ranks(self) -> None:
        """Update model_rank feature for all checkpoints under this key."""
        artifacts = ln.Artifact.filter(**self._get_key_filter()).to_list()
        scored = []
        for af in artifacts:
            vals = af.features.get_values()
            if "score" in vals:
                scored.append((vals["score"], vals.get("model_rank"), af))
        scored.sort(key=lambda x: x[0], reverse=(self.mode == "max"))
        for rank, (_, old_rank, af) in enumerate(scored):
            if old_rank is not None:
                af.features.remove_values("model_rank", value=old_rank)
            af.features.add_values({"model_rank": rank})


class SaveConfigCallback(_SaveConfigCallback):
    """SaveConfigCallback that also saves config to the instance.

    Use with LightningCLI to save the resolved configuration file alongside checkpoints.

    Example::

        from lightning.pytorch.cli import LightningCLI
        from lamindb.integrations import lightning as ll

        cli = LightningCLI(
            MyModel,
            MyDataModule,
            save_config_callback=ll.SaveConfigCallback,
        )
    """

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Save resolved configuration file alongside checkpoints."""
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
                self._save_config(trainer, config_path)

            if trainer.is_global_zero:
                self.save_config(trainer, pl_module, stage)
                self.already_saved = True
            self.already_saved = trainer.strategy.broadcast(self.already_saved)

    def _save_config(self, trainer: pl.Trainer, config_path: Path) -> None:
        """Save config under same key prefix as checkpoints."""
        checkpoint_cb = self._get_checkpoint_callback(trainer)
        if checkpoint_cb is None:
            return

        lightning_cli_config_af = ln.Artifact(
            config_path,
            key=f"{checkpoint_cb.dirpath.rstrip('/')}/{self.config_filename}",
            description="Lightning CLI config",
        )
        lightning_cli_config_af.save()

        # Link as input to current run
        if ln.context.run:
            ln.context.run.input_artifacts.add(lightning_cli_config_af)

    def _get_checkpoint_callback(self, trainer: pl.Trainer) -> Checkpoint | None:
        """Find LaminCheckpoint callback if present."""
        for cb in trainer.callbacks:
            if isinstance(cb, Checkpoint):
                return cb
        return None


# backwards compatibility
# We keep the full class around because it's short and it's cumbersome to write
# full backwards compatibility code because of the rather different interfaces and behavior
class Callback(pl.Callback):
    """Saves checkpoints to LaminDB after each training epoch.

    .. deprecated::
        Use :class:`Checkpoint` instead for new code.

    Args:
        path: A local path to the checkpoint.
        key: The `key` for the checkpoint artifact.
        features: Features to annotate the checkpoint.
    """

    def __init__(
        self,
        path: str | Path,
        key: str,
        features: dict[str, Any] | None = None,
    ):
        warnings.warn(
            "ll.Callback is deprecated, use ll.Checkpoint instead",
            DeprecationWarning,
            stacklevel=2,
        )
        self.path = Path(path)
        self.key = key
        self.features = features or {}

    def on_train_start(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule
    ) -> None:
        """Validates that features exist for all specified params."""
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

    def on_train_epoch_end(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule
    ) -> None:
        """Saves model checkpoint at the end of each epoch."""
        trainer.save_checkpoint(self.path)
        artifact = ln.Artifact(self.path, key=self.key, kind="model").save()

        feature_values = dict(self.features)
        for name in self.features:
            if hasattr(trainer, name):
                feature_values[name] = getattr(trainer, name)
            elif name in trainer.callback_metrics:
                metric = trainer.callback_metrics[name]
                feature_values[name] = (
                    metric.item() if hasattr(metric, "item") else float(metric)
                )

        if feature_values:
            artifact.features.add_values(feature_values)


__all__ = ["Checkpoint", "SaveConfigCallback", "save_lightning_features"]
