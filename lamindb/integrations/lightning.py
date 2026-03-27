"""PyTorch Lightning integration for LaminDB.

The public surface has two layers:

- :class:`Checkpoint` is the concrete LaminDB implementation that persists
    checkpoint, config, and ``hparams.yaml`` files as :class:`lamindb.Artifact`
    records and annotates them with Lamin features.
- :class:`ArtifactPublishingModelCheckpoint` is the generic extension layer adding
    checkpoint artifact lifecycle hooks without implementing Lamin persistence
    details yet. Mostly done for clarity.

External integrations can either subclass :class:`Checkpoint` directly or attach
an :class:`ArtifactObserver` to react to saved and removed artifacts.

.. autoclass:: ArtifactPublishingModelCheckpoint
.. autoclass:: Checkpoint
.. autoclass:: SaveConfigCallback
.. autoclass:: ArtifactSavedEvent
.. autoclass:: ArtifactRemovedEvent
.. autofunction:: save_lightning_features
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Literal, Protocol

import lightning.pytorch as pl
from lightning.fabric.utilities.cloud_io import get_filesystem
from lightning.pytorch.callbacks.model_checkpoint import ModelCheckpoint
from lightning.pytorch.cli import SaveConfigCallback as _SaveConfigCallback
from lamindb.models.artifact import track_run_input

import lamindb as ln

if TYPE_CHECKING:
    from datetime import timedelta

    from lightning.fabric.utilities.types import _PATH


_RUN_AUTO_FEATURES: Final = frozenset(
    {
        "logger_name",
        "logger_version",
        "max_epochs",
        "max_steps",
        "precision",
        "accumulate_grad_batches",
        "gradient_clip_val",
        "monitor",
        "save_weights_only",
        "mode",
    }
)
_ARTIFACT_AUTO_FEATURES: Final = frozenset({"is_best_model", "score", "model_rank"})
_SUPPORTED_AUTO_FEATURES: Final = _RUN_AUTO_FEATURES | _ARTIFACT_AUTO_FEATURES
ArtifactKind = Literal["checkpoint", "config", "hparams"]


@dataclass(frozen=True)
class ArtifactEvent:
    """Common metadata emitted when a checkpoint-related artifact changes.

    The event records the logical artifact key, the local path Lightning wrote,
    and the trainer that triggered the lifecycle event.
    """

    kind: ArtifactKind
    key: str
    local_path: Path
    trainer: pl.Trainer


@dataclass(frozen=True)
class ArtifactSavedEvent(ArtifactEvent):
    """Metadata emitted after a checkpoint-related artifact has been persisted.

    ``artifact`` is intentionally typed generically so downstream integrations can
    expose their own persisted object while still using the common lifecycle API.
    ``storage_uri`` is the stable hand-off value for registries such as ClearML.
    """

    artifact: Any
    storage_uri: str


@dataclass(frozen=True)
class ArtifactRemovedEvent(ArtifactEvent):
    """Metadata emitted after a local checkpoint file has been removed.

    Removal currently applies to checkpoint files. Config and hparams artifacts are
    save-only in the current Lightning integration.
    """

    artifact: Any | None = None
    storage_uri: str | None = None


class ArtifactObserver(Protocol):
    """Observer notified about checkpoint artifact lifecycle events.

    This is the preferred composition hook for downstream integrations that need
    to register checkpoints elsewhere after Lamin persistence completes.
    """

    def on_artifact_saved(self, event: ArtifactSavedEvent) -> None: ...

    def on_artifact_removed(self, event: ArtifactRemovedEvent) -> None: ...


class ArtifactPublisher(Protocol):
    """Persistence backend for checkpoint-related artifacts.

    :class:`ArtifactPublishingModelCheckpoint` manages the artifact lifecycle,
    while publishers encapsulate backend-specific save behavior and storage URI
    resolution.
    """

    def create_artifact(
        self,
        local_path: Path | str,
        *,
        key: str,
        description: str,
        kind: str | None = None,
        add_as_input_to_run: bool = False,
        skip_hash_lookup: bool = False,
    ) -> Any: ...

    def storage_uri(self, artifact: Any) -> str: ...


class LaminArtifactPublisher:
    """Persist checkpoint-related artifacts into LaminDB.

    This service is intentionally separate from :class:`Checkpoint` so that the
    checkpoint callback can focus on Lightning behavior and feature handling while
    persistence details remain replaceable.
    """

    def create_artifact(
        self,
        local_path: Path | str,
        *,
        key: str,
        description: str,
        kind: str | None = None,
        add_as_input_to_run: bool = False,
        skip_hash_lookup: bool = False,
    ) -> ln.Artifact:
        artifact_kwargs: dict[str, Any] = {"key": key, "description": description}
        if kind is not None:
            artifact_kwargs["kind"] = kind
        if add_as_input_to_run:
            artifact_kwargs["run"] = False
        if skip_hash_lookup:
            artifact_kwargs["skip_hash_lookup"] = True
        artifact = ln.Artifact(local_path, **artifact_kwargs)
        artifact.save()
        if add_as_input_to_run:
            track_run_input(artifact, is_run_input=True)
        return artifact

    def storage_uri(self, artifact: ln.Artifact) -> str:
        return str(artifact.path)


def save_lightning_features() -> None:
    """Save features to auto-track lightning parameters & metrics.

    Creates the following features under the `lamindb.lightning` feature type if they do not already exist:

    - `is_best_model` (bool): Whether this checkpoint is the best model.
    - `score` (float): The monitored metric score.
    - `model_rank` (int): Rank among all checkpoints (0 = best).
    - `logger_name` (str): Name from the first Lightning logger.
    - `logger_version` (str): Version from the first Lightning logger.
    - `max_epochs` (int): Maximum number of epochs.
    - `max_steps` (int): Maximum number of training steps.
    - `precision` (str): Training precision (e.g., "32", "16-mixed", "bf16").
    - `accumulate_grad_batches` (int): Number of batches to accumulate gradients over.
    - `gradient_clip_val` (float): Gradient clipping value.
    - `monitor` (str): Metric name being monitored.
    - `save_weights_only` (bool): Whether only model weights are saved.
    - `mode` (str): Optimization mode ("min" or "max").

    Args:
        None.

    Example:

        Save the features to the database::

            from lamindb.integrations import lightning as ll

            ll.save_lightning_features()
    """
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
    ln.Feature(name="max_epochs", dtype=int, type=lightning_feature_type).save()
    ln.Feature(name="max_steps", dtype=int, type=lightning_feature_type).save()
    ln.Feature(name="precision", dtype=str, type=lightning_feature_type).save()
    ln.Feature(
        name="accumulate_grad_batches", dtype=int, type=lightning_feature_type
    ).save()
    ln.Feature(
        name="gradient_clip_val", dtype=float, type=lightning_feature_type
    ).save()
    ln.Feature(name="monitor", dtype=str, type=lightning_feature_type).save()
    ln.Feature(name="save_weights_only", dtype=bool, type=lightning_feature_type).save()
    ln.Feature(name="mode", dtype=str, type=lightning_feature_type).save()


class ArtifactPublishingModelCheckpoint(ModelCheckpoint):
    """ModelCheckpoint with observable artifact lifecycle hooks.

    This layer captures artifact kinds, observer registration, saved/removed
        events, latest artifact tracking, and key compatibility hooks. Concrete
        subclasses remain responsible for how artifacts are persisted.

        Subclasses are expected to implement:

        - :meth:`resolve_artifact_key` to map local files to logical artifact keys
        - :meth:`resolve_artifact_storage_uri` to expose a stable backend URI
        - :meth:`save_checkpoint_artifact`, :meth:`save_config_artifact`, and
            :meth:`save_hparams_artifact` to persist files

        :class:`SaveConfigCallback` only depends on this base class, which means a
        custom checkpoint callback can participate in config saving without inheriting
        from Lamin's concrete :class:`Checkpoint`.
    """

    def __init__(
        self,
        *args: Any,
        artifact_observers: list[ArtifactObserver] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self._artifact_observers: list[ArtifactObserver] = list(
            artifact_observers or []
        )
        self._latest_artifacts: dict[ArtifactKind, Any | None] = {
            "checkpoint": None,
            "config": None,
            "hparams": None,
        }
        self._last_artifact_event: ArtifactSavedEvent | ArtifactRemovedEvent | None = (
            None
        )

    @property
    def last_checkpoint_artifact(self) -> Any | None:
        """The most recently saved checkpoint artifact handle."""
        return self._latest_artifacts["checkpoint"]

    @property
    def last_config_artifact(self) -> Any | None:
        """The most recently saved config artifact handle."""
        return self._latest_artifacts["config"]

    @property
    def last_hparams_artifact(self) -> Any | None:
        """The most recently saved hparams artifact handle."""
        return self._latest_artifacts["hparams"]

    @property
    def last_artifact_event(self) -> ArtifactSavedEvent | ArtifactRemovedEvent | None:
        """The last artifact lifecycle event emitted by this callback."""
        return self._last_artifact_event

    def get_last_artifact(self, kind: ArtifactKind) -> Any | None:
        """Return the most recently saved artifact for a given artifact kind."""
        return self._latest_artifacts[kind]

    def add_artifact_observer(self, observer: ArtifactObserver) -> None:
        """Register an observer notified about artifact lifecycle events."""
        self._artifact_observers.append(observer)

    def remove_artifact_observer(self, observer: ArtifactObserver) -> None:
        """Unregister a previously added artifact observer."""
        self._artifact_observers.remove(observer)

    def resolve_artifact_storage_uri(self, artifact: Any) -> str:
        """Resolve the physical location for a persisted artifact."""
        raise NotImplementedError

    def resolve_artifact_key(
        self,
        trainer: pl.Trainer,
        filepath: Path | str,
        kind: ArtifactKind,
    ) -> str:
        """Return the logical artifact key for a checkpoint-related file."""
        raise NotImplementedError

    def _get_artifact_key(
        self, trainer: pl.Trainer, filepath: Path | str, is_checkpoint: bool = True
    ) -> str:
        """Backward-compatible wrapper around :meth:`resolve_artifact_key`."""
        kind: ArtifactKind = "checkpoint" if is_checkpoint else "config"
        return self.resolve_artifact_key(trainer=trainer, filepath=filepath, kind=kind)

    def _artifact_key_for_kind(
        self, trainer: pl.Trainer, filepath: Path | str, kind: ArtifactKind
    ) -> str:
        """Return the artifact key while preserving legacy checkpoint overrides."""
        if kind == "checkpoint":
            return self._get_artifact_key(
                trainer=trainer, filepath=filepath, is_checkpoint=True
            )
        if kind == "config":
            return self._get_artifact_key(
                trainer=trainer, filepath=filepath, is_checkpoint=False
            )
        return self.resolve_artifact_key(trainer=trainer, filepath=filepath, kind=kind)

    def _notify_artifact_saved(
        self,
        trainer: pl.Trainer,
        *,
        kind: ArtifactKind,
        artifact: Any,
        local_path: Path | str,
    ) -> ArtifactSavedEvent:
        event = ArtifactSavedEvent(
            kind=kind,
            key=artifact.key,
            local_path=Path(local_path),
            trainer=trainer,
            artifact=artifact,
            storage_uri=self.resolve_artifact_storage_uri(artifact),
        )
        self._latest_artifacts[kind] = artifact
        self._last_artifact_event = event
        self.on_artifact_saved(event)
        self._notify_artifact_observers("on_artifact_saved", event)
        return event

    def _notify_artifact_removed(
        self,
        trainer: pl.Trainer,
        *,
        kind: ArtifactKind,
        key: str,
        local_path: Path | str,
        artifact: Any | None,
    ) -> ArtifactRemovedEvent:
        storage_uri = None
        if artifact is not None:
            storage_uri = self.resolve_artifact_storage_uri(artifact)
        event = ArtifactRemovedEvent(
            kind=kind,
            key=key,
            local_path=Path(local_path),
            trainer=trainer,
            artifact=artifact,
            storage_uri=storage_uri,
        )
        self._last_artifact_event = event
        self.on_artifact_removed(event)
        self._notify_artifact_observers("on_artifact_removed", event)
        return event

    def _notify_artifact_observers(
        self,
        method_name: str,
        event: ArtifactSavedEvent | ArtifactRemovedEvent,
    ) -> None:
        for observer in tuple(self._artifact_observers):
            method = getattr(observer, method_name, None)
            if callable(method):
                method(event)

    def on_artifact_saved(self, event: ArtifactSavedEvent) -> None:
        """Hook for subclasses after an artifact has been saved."""
        del event

    def on_artifact_removed(self, event: ArtifactRemovedEvent) -> None:
        """Hook for subclasses after a checkpoint file has been removed."""
        del event

    def save_checkpoint_artifact(
        self,
        trainer: pl.Trainer,
        filepath: Path | str,
        *,
        feature_values: dict[str, Any] | None = None,
    ) -> Any:
        """Persist a checkpoint artifact and emit the corresponding event."""
        del trainer, filepath, feature_values
        raise NotImplementedError

    def save_config_artifact(self, trainer: pl.Trainer, config_path: Path | str) -> Any:
        """Persist a config artifact and emit the corresponding event."""
        del trainer, config_path
        raise NotImplementedError

    def save_hparams_artifact(
        self, trainer: pl.Trainer, hparams_path: Path | str
    ) -> Any | None:
        """Persist an hparams artifact and emit the corresponding event."""
        del trainer, hparams_path
        raise NotImplementedError


class Checkpoint(ArtifactPublishingModelCheckpoint):
    """A `ModelCheckpoint` that annotates `pytorch` `lightning` checkpoints.

    Extends `lightning`'s `ModelCheckpoint` with artifact creation & feature annotation.
    Each checkpoint is a separate artifact whose key is derived from either the
    explicit ``dirpath`` or the trainer's logger configuration.

    When ``dirpath`` is omitted (recommended), Lightning decides where to store
    checkpoints locally (typically ``lightning_logs/version_N/checkpoints/``)
    and the artifact key is derived from the logger's ``save_dir``, ``name``,
    and ``version``.  When ``dirpath`` is provided, it is used directly as the
    key prefix.

    Key derivation rules:

    - ``dirpath`` provided → ``{dirpath}/{filename}``
    - ``dirpath`` omitted, logger present → ``{save_dir_basename}/{name}/{version}/checkpoints/{filename}``
    - ``dirpath`` omitted, no logger → ``checkpoints/{filename}``

    If available in the database through `save_lightning_features()`, the following `lamindb.lightning` features are automatically tracked:
    `is_best_model`, `score`, `model_rank`, `logger_name`, `logger_version`,`max_epochs`, `max_steps`,
    `precision`, `accumulate_grad_batches`, `gradient_clip_val`, `monitor`, `save_weights_only`, `mode`.

    Additionally, model hyperparameters (from `pl_module.hparams`) and datamodule hyperparameters
    (from `trainer.datamodule.hparams`) are captured if corresponding features exist.

    This is the concrete LaminDB implementation built on top of
    :class:`ArtifactPublishingModelCheckpoint`. Use it when you want LaminDB to be
    the persistence layer. For secondary systems such as ClearML, prefer attaching
    an :class:`ArtifactObserver` or subclassing :class:`Checkpoint` and reacting in
    :meth:`on_artifact_saved`.

    Args:
        dirpath: Directory for checkpoints.  When provided, also used as the
            artifact key prefix.  When omitted (recommended), Lightning picks
            the local directory and the key prefix is derived from the logger.
        features: Features to annotate runs and artifacts.
            Use "run" key for run-level features (static metadata).
            Use "artifact" key for artifact-level features (values can be static or None for auto-population from trainer metrics/attributes).
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
        overwrite_versions: Whether to overwrite existing checkpoints.
        artifact_observers: Optional observer objects notified when checkpoint,
            config, or hparams artifacts are saved or when checkpoint files are
            removed locally. Observers follow :class:`ArtifactObserver` and
            receive :class:`ArtifactSavedEvent` and :class:`ArtifactRemovedEvent`.

    Examples:

        Let Lightning decide where to store checkpoints (recommended)::

            import lightning as pl
            from lightning.pytorch.loggers import CSVLogger
            from lamindb.integrations import lightning as ll

            ll.save_lightning_features()

            callback = ll.Checkpoint(monitor="val_loss", save_top_k=3)
            logger = CSVLogger(save_dir="logs")

            trainer = pl.Trainer(callbacks=[callback], logger=logger)
            trainer.fit(model, dataloader)

            # Query checkpoints — key prefix is derived from the logger
            # e.g. "logs/lightning_logs/version_0/checkpoints/"
            ln.Artifact.filter(key__startswith=callback.checkpoint_key_prefix)

        Explicit ``dirpath`` for full control over the artifact key prefix::

            callback = ll.Checkpoint(
                dirpath="deployments/my_model/",
                monitor="val_loss",
                save_top_k=3,
            )

            trainer = pl.Trainer(callbacks=[callback])
            trainer.fit(model, dataloader)

            # Query checkpoints
            ln.Artifact.filter(key__startswith=callback.checkpoint_key_prefix)

        Using the CLI::

            # config.yaml
            trainer:
              callbacks:
                - class_path: lamindb.integrations.lightning.Checkpoint
                  init_args:
                    monitor: val_loss
                    save_top_k: 3

            # Run with:
            # python main.py fit --config config.yaml
    """

    def __init__(
        self,
        dirpath: _PATH | None = None,
        *,
        features: dict[Literal["run", "artifact"], dict[str, Any]] | None = None,
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
        overwrite_versions: bool = False,
        artifact_observers: list[ArtifactObserver] | None = None,
    ) -> None:
        self._original_dirpath = dirpath
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
            artifact_observers=artifact_observers,
        )
        self.features = features or {}
        if invalid_feature_keys := set(self.features) - {"run", "artifact"}:  # type: ignore
            raise ValueError(
                f"Invalid feature keys: {invalid_feature_keys}. Use 'run' and/or 'artifact'."
            )
        self._run_features = self.features.get("run", {})
        self._artifact_features = self.features.get("artifact", {})
        self._auto_features: dict[str, ln.Feature] = {}
        self._hparam_features_available: set[str] = set()
        self._run_features_added = False
        self._hparams_yaml_saved = False
        self.overwrite_versions = overwrite_versions
        self._checkpoint_key_prefix: str = ""
        self._artifact_publisher: ArtifactPublisher = LaminArtifactPublisher()

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Validate user features and detect available auto-features."""
        super().setup(trainer, pl_module, stage)
        self._checkpoint_key_prefix = self._compute_checkpoint_key_prefix(trainer)

        if trainer.is_global_zero:
            # Validate user-specified features exist
            all_feature_names = set(self._run_features) | set(self._artifact_features)
            existing_feature_names = set(
                ln.Feature.filter(name__in=all_feature_names).values_list(
                    "name", flat=True
                )
            )
            missing = [
                name for name in all_feature_names if name not in existing_feature_names
            ]
            if missing:
                s = "s" if len(missing) > 1 else ""
                raise ValueError(
                    f"Feature{s} {', '.join(missing)} missing. "
                    f"Create {'them' if len(missing) > 1 else 'it'} first."
                )

            if ln.context.run and self._run_features:
                ln.context.run.features.add_values(self._run_features)

            # Auto-features are opt-in and scoped to the lightning feature type.
            # If `save_lightning_features()` was never called, skip auto-tracking.
            lightning_feature_type = ln.Feature.filter(
                name="lamindb.lightning", is_type=True
            ).one_or_none()
            self._auto_features.clear()
            if lightning_feature_type is not None:
                auto_features = list(
                    ln.Feature.filter(
                        name__in=_SUPPORTED_AUTO_FEATURES,
                        type=lightning_feature_type,
                    )
                )
                self._auto_features = {f.name: f for f in auto_features}
            hparam_names: set[str] = set()
            if hasattr(pl_module, "hparams") and pl_module.hparams:
                hparam_names.update(pl_module.hparams.keys())
            if (
                trainer.datamodule is not None
                and hasattr(trainer.datamodule, "hparams")
                and trainer.datamodule.hparams
            ):
                hparam_names.update(trainer.datamodule.hparams.keys())
            self._hparam_features_available = set(
                ln.Feature.filter(name__in=hparam_names).values_list("name", flat=True)
            )

    @staticmethod
    def _compute_logger_key_prefix(trainer: pl.Trainer) -> str:
        """Derive a key prefix from the trainer's first logger."""
        if trainer.loggers[0].save_dir is not None:
            save_dir = trainer.loggers[0].save_dir
        else:
            save_dir = trainer.default_root_dir
        name = trainer.loggers[0].name
        version = trainer.loggers[0].version
        version = version if isinstance(version, str) else f"version_{version}"
        return f"{Path(save_dir).name}/{str(name).rstrip('/')}/{version.rstrip('/')}"

    def _compute_checkpoint_key_prefix(self, trainer: pl.Trainer) -> str:
        """Compute the artifact key prefix for checkpoints."""
        if self._original_dirpath is not None:
            return str(self._original_dirpath).rstrip("/")
        if len(trainer.loggers) > 0:
            return f"{self._compute_logger_key_prefix(trainer)}/checkpoints"
        return "checkpoints"

    @property
    def checkpoint_key_prefix(self) -> str:
        """The artifact key prefix used for checkpoint artifacts.

        Available after ``setup()`` has been called, for example once
        ``trainer.fit()`` has started.
        """
        return self._checkpoint_key_prefix

    def resolve_artifact_storage_uri(self, artifact: ln.Artifact) -> str:
        """Resolve the physical artifact location for downstream registries.

        This is the stable abstraction external packages should use instead of
        reconstructing storage locations from Lamin internals.
        """
        return self._artifact_publisher.storage_uri(artifact)

    def resolve_artifact_key(
        self,
        trainer: pl.Trainer,
        filepath: Path | str,
        kind: ArtifactKind,
    ) -> str:
        """Return the Lamin artifact key for a checkpoint-related file."""
        if kind in {"checkpoint", "hparams"}:
            prefix = self._checkpoint_key_prefix
        elif len(trainer.loggers) > 0:
            prefix = self._compute_logger_key_prefix(trainer)
        elif self._original_dirpath is not None:
            prefix = str(self._original_dirpath).rstrip("/")
        else:
            prefix = ""
        if prefix:
            return f"{prefix}/{Path(filepath).name}"
        return Path(filepath).name

    def _get_key_filter(self) -> dict[str, str]:
        """Return filter kwargs for querying artifacts from this callback."""
        return {"key__startswith": self._checkpoint_key_prefix + "/"}

    def _create_lamin_artifact(
        self,
        local_path: Path | str,
        *,
        key: str,
        description: str,
        kind: str | None = None,
        add_as_input_to_run: bool = False,
        skip_hash_lookup: bool = False,
    ) -> ln.Artifact:
        return self._artifact_publisher.create_artifact(
            local_path,
            key=key,
            description=description,
            kind=kind,
            add_as_input_to_run=add_as_input_to_run,
            skip_hash_lookup=skip_hash_lookup,
        )

    def save_checkpoint_artifact(
        self,
        trainer: pl.Trainer,
        filepath: Path | str,
        *,
        feature_values: dict[str | ln.Feature, Any] | None = None,
    ) -> ln.Artifact:
        """Save a checkpoint artifact to Lamin and emit the corresponding event.

        This is the main persistence hook used by :meth:`_save_checkpoint`. It is a
        useful override point for subclasses that want to augment Lamin persistence
        while keeping the generic lifecycle behavior from the base class.
        """
        key = self._artifact_key_for_kind(
            trainer=trainer, filepath=filepath, kind="checkpoint"
        )
        existing_artifact = ln.Artifact.filter(key=key).one_or_none()
        if existing_artifact is not None:
            if not self.overwrite_versions:
                raise ValueError(
                    f"An artifact with key '{key}' already exists. Set overwrite_versions=True to replace it."
                )
            existing_artifact.delete(permanent=True, storage=True)
        artifact = self._create_lamin_artifact(
            filepath,
            key=key,
            description="model checkpoint",
            kind="model",
            skip_hash_lookup=True,
        )
        if feature_values:
            artifact.features.add_values(feature_values)
        self._notify_artifact_saved(
            trainer,
            kind="checkpoint",
            artifact=artifact,
            local_path=filepath,
        )
        return artifact

    def save_config_artifact(
        self, trainer: pl.Trainer, config_path: Path | str
    ) -> ln.Artifact:
        """Save a Lightning CLI config artifact and emit the corresponding event.

        Config artifacts are routed through the same lifecycle surface as
        checkpoints so observers and subclasses see a unified event stream.
        """
        key = self._artifact_key_for_kind(
            trainer=trainer, filepath=config_path, kind="config"
        )
        artifact = self._create_lamin_artifact(
            config_path,
            key=key,
            description="Lightning CLI config",
            kind="config",
            add_as_input_to_run=True,
            skip_hash_lookup=True,
        )
        self._notify_artifact_saved(
            trainer,
            kind="config",
            artifact=artifact,
            local_path=config_path,
        )
        return artifact

    def save_hparams_artifact(
        self, trainer: pl.Trainer, hparams_path: Path | str
    ) -> ln.Artifact | None:
        """Save Lightning's auto-generated hparams file and emit the event.

        Returns ``None`` if Lightning did not generate ``hparams.yaml`` for the
        current run.
        """
        if not Path(hparams_path).exists():
            return None

        key = self._artifact_key_for_kind(
            trainer=trainer, filepath=hparams_path, kind="hparams"
        )
        artifact = self._create_lamin_artifact(
            hparams_path,
            key=key,
            description="Lightning run hyperparameters",
            kind="config",
            skip_hash_lookup=True,
        )
        self._notify_artifact_saved(
            trainer,
            kind="hparams",
            artifact=artifact,
            local_path=hparams_path,
        )
        return artifact

    def _save_hparams_yaml(self, trainer: pl.Trainer) -> None:
        """Persist Lightning's auto-generated hparams file once per run."""
        if self._hparams_yaml_saved:
            return

        log_dir = trainer.log_dir
        if not log_dir:
            return

        hparams_path = Path(log_dir) / "hparams.yaml"
        if not hparams_path.exists():
            return

        if self.save_hparams_artifact(trainer, hparams_path) is not None:
            self._hparams_yaml_saved = True

    def _save_run_features(self, trainer: pl.Trainer) -> None:
        """Save run-level features once per run after the first checkpoint."""
        if not (ln.context.run and not self._run_features_added):
            return

        run_features: dict[str | ln.Feature, Any] = {}
        if (feature := self._get_lightning_feature("logger_name")) and trainer.loggers:
            run_features[feature] = trainer.loggers[0].name
        if (feature := self._get_lightning_feature("logger_version")) and trainer.loggers:
            version = trainer.loggers[0].version
            run_features[feature] = (
                version if isinstance(version, str) else f"version_{version}"
            )
        trainer_config_values = {
            "max_epochs": (trainer.max_epochs, None, True),
            "max_steps": (trainer.max_steps, None, True),
            "precision": (trainer.precision, str, False),
            "accumulate_grad_batches": (
                trainer.accumulate_grad_batches,
                None,
                False,
            ),
            "gradient_clip_val": (trainer.gradient_clip_val, float, True),
            "monitor": (self.monitor, None, True),
            "save_weights_only": (self.save_weights_only, None, False),
            "mode": (self.mode, None, False),
        }

        for key, (value, transform, skip_none) in trainer_config_values.items():
            if (feature := self._get_lightning_feature(key)) and not (
                skip_none and value is None
            ):
                run_features[feature] = transform(value) if transform else value

        if (
            hasattr(trainer.lightning_module, "hparams")
            and trainer.lightning_module.hparams
        ):
            for name, value in trainer.lightning_module.hparams.items():
                if name in self._hparam_features_available:
                    run_features[name] = value
        if (
            trainer.datamodule is not None
            and hasattr(trainer.datamodule, "hparams")
            and trainer.datamodule.hparams
        ):
            for name, value in trainer.datamodule.hparams.items():
                if name in self._hparam_features_available:
                    run_features[name] = value
        if run_features:
            ln.context.run.features.add_values(run_features)
        self._run_features_added = True

    def _collect_checkpoint_feature_values(
        self, trainer: pl.Trainer, filepath: Path | str
    ) -> dict[str | ln.Feature, Any]:
        """Collect checkpoint-level feature values before saving the artifact."""
        feature_values: dict[str | ln.Feature, Any] = {}
        is_best = self.best_model_path == str(filepath)
        if feature := self._get_lightning_feature("is_best_model"):
            if is_best:
                self._clear_best_model_flags()
            feature_values[feature] = is_best

        import torch

        if (feature := self._get_lightning_feature("score")) and self.current_score is not None:
            score = self.current_score
            if torch.is_tensor(score):
                score = score.item()
            feature_values[feature] = float(score)

        for name, value in self._artifact_features.items():
            if value is not None:
                feature_values[name] = value
            elif hasattr(trainer, name):
                feature_values[name] = getattr(trainer, name)
            elif name in trainer.callback_metrics:
                metric = trainer.callback_metrics[name]
                feature_values[name] = (
                    metric.item() if hasattr(metric, "item") else float(metric)
                )
        return feature_values

    def _save_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Save checkpoint to the instance."""
        super()._save_checkpoint(trainer, filepath)

        if trainer.is_global_zero:
            self._save_hparams_yaml(trainer)

            self._save_run_features(trainer)
            feature_values = self._collect_checkpoint_feature_values(trainer, filepath)
            self.save_checkpoint_artifact(
                trainer, filepath, feature_values=feature_values
            )

            if self._get_lightning_feature("model_rank"):
                self._update_model_ranks()

    def _remove_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Remove the local checkpoint file and emit a removal event."""
        artifact: ln.Artifact | None = None
        key = self._artifact_key_for_kind(
            trainer=trainer, filepath=filepath, kind="checkpoint"
        )
        if trainer.is_global_zero:
            artifact = ln.Artifact.filter(key=key).one_or_none()
        super()._remove_checkpoint(trainer, filepath)
        if trainer.is_global_zero:
            self._notify_artifact_removed(
                trainer,
                kind="checkpoint",
                key=key,
                local_path=filepath,
                artifact=artifact,
            )

    def _get_lightning_feature(self, name: str) -> ln.Feature | None:
        """Return the cached typed Feature for *name*, or None."""
        return self._auto_features.get(name)

    def _clear_best_model_flags(self) -> None:
        """Set is_best_model=False on previous best checkpoints."""
        is_best_feature = self._get_lightning_feature("is_best_model")
        if is_best_feature is None:
            return
        feature_rows = self._get_artifact_feature_rows({"is_best_model"})
        best_artifact_ids = [
            artifact_id
            for artifact_id, values in feature_rows.items()
            if values.get("is_best_model") is True
        ]
        if not best_artifact_ids:
            return
        artifacts_by_id = {
            artifact.id: artifact
            for artifact in ln.Artifact.filter(id__in=best_artifact_ids)
        }
        for artifact_id in best_artifact_ids:
            if artifact_id not in artifacts_by_id:
                continue
            artifact = artifacts_by_id[artifact_id]
            artifact.features.remove_values(is_best_feature, value=True)
            artifact.features.add_values({is_best_feature: False})

    def _update_model_ranks(self) -> None:
        """Update model_rank feature for all checkpoints under this key."""
        model_rank_feature = self._get_lightning_feature("model_rank")
        if model_rank_feature is None:
            return
        feature_rows = self._get_artifact_feature_rows({"score", "model_rank"})
        scored = []
        for artifact_id, values in feature_rows.items():
            if "score" in values:
                scored.append((values["score"], values.get("model_rank"), artifact_id))
        scored.sort(key=lambda x: x[0], reverse=(self.mode == "max"))

        artifact_ids = [artifact_id for _, _, artifact_id in scored]
        artifacts_by_id = {
            artifact.id: artifact
            for artifact in ln.Artifact.filter(id__in=artifact_ids)
        }
        for rank, (_, old_rank, artifact_id) in enumerate(scored):
            if artifact_id not in artifacts_by_id:
                continue
            af = artifacts_by_id[artifact_id]
            if old_rank is not None:
                af.features.remove_values(model_rank_feature, value=old_rank)
            af.features.add_values({model_rank_feature: rank})

    def _get_artifact_feature_rows(
        self, feature_names: set[str]
    ) -> dict[int, dict[str, Any]]:
        # Filter by feature id when possible to avoid cross-type name collisions
        feature_ids = [
            feature.id
            for name in feature_names
            if (feature := self._get_lightning_feature(name))
        ]
        key_prefix = self._get_key_filter()["key__startswith"]
        if feature_ids:
            rows = ln.models.ArtifactJsonValue.filter(
                artifact__key__startswith=key_prefix,
                jsonvalue__feature_id__in=feature_ids,
            ).values_list(
                "artifact_id", "jsonvalue__feature__name", "jsonvalue__value"
            )
        else:
            rows = ln.models.ArtifactJsonValue.filter(
                artifact__key__startswith=key_prefix,
                jsonvalue__feature__name__in=feature_names,
            ).values_list(
                "artifact_id", "jsonvalue__feature__name", "jsonvalue__value"
            )
        result: dict[int, dict[str, Any]] = {}
        for artifact_id, feature_name, value in rows:
            if artifact_id not in result:
                result[artifact_id] = {}
            result[artifact_id][feature_name] = value
        return result


class SaveConfigCallback(_SaveConfigCallback):
    """SaveConfigCallback that also saves config to the instance.

    Use with LightningCLI to save the resolved configuration file alongside checkpoints.

    The local config file follows Lightning's ``trainer.log_dir`` behavior.
    Setting ``dirpath`` on :class:`Checkpoint` affects where checkpoint files are
    saved, but normally does not affect the config file location.

    This callback looks for any :class:`ArtifactPublishingModelCheckpoint`, not just
    Lamin's concrete :class:`Checkpoint`. That keeps the config-save path aligned
    with custom subclasses built on the generic artifact-publishing base.

    Lamin stores config artifacts with the following key rules:

    - logger present -> ``{save_dir_basename}/{name}/{version}/config.yaml``
    - no logger, ``Checkpoint.dirpath`` set -> ``{dirpath}/config.yaml``
    - no logger, no ``dirpath`` -> ``config.yaml``

    The no-logger + explicit-``dirpath`` case is a Lamin-specific exception to avoid
    storing ``config.yaml`` at the storage root when checkpoints already have an
    explicit namespace. The local config file still remains under ``trainer.log_dir``.

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
        """Persist the resolved config through the active artifact checkpoint.

        If no artifact-publishing checkpoint callback is registered, this becomes a
        no-op and only Lightning's local config file is written.
        """
        checkpoint_cb = self._get_artifact_checkpoint_callback(trainer)
        if checkpoint_cb is None:
            return

        checkpoint_cb.save_config_artifact(trainer, config_path)

    def _get_artifact_checkpoint_callback(
        self, trainer: pl.Trainer
    ) -> ArtifactPublishingModelCheckpoint | None:
        """Find the artifact-publishing checkpoint callback if present."""
        for cb in trainer.callbacks:
            if isinstance(cb, ArtifactPublishingModelCheckpoint):
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


__all__ = [
    "ArtifactObserver",
    "ArtifactEvent",
    "ArtifactPublisher",
    "ArtifactPublishingModelCheckpoint",
    "ArtifactRemovedEvent",
    "ArtifactSavedEvent",
    "Checkpoint",
    "LaminArtifactPublisher",
    "SaveConfigCallback",
    "save_lightning_features",
]
