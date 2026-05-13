"""PyTorch Lightning integration for LaminDB.

The public API has two layers:

- :class:`Checkpoint` is the concrete LaminDB implementation that persists checkpoint, config, and `hparams.yaml` files as :class:`~lamindb.Artifact` objects and annotates them with :class:`~lamindb.Feature` objects.
- :class:`ArtifactPublishingModelCheckpoint` is the generic extension layer adding checkpoint artifact lifecycle hooks without implementing Lamin persistence details yet.

External integrations can either subclass :class:`Checkpoint` directly or attach
an :class:`ArtifactObserver` to react to saved and removed artifacts.

Here is a guide: :doc:`lightning`.

Main API
--------

.. autoclass:: Checkpoint
.. autofunction:: save_lightning_features

Auxiliary classes
-----------------

.. autoclass:: ArtifactPublishingModelCheckpoint
.. autoclass:: SaveConfigCallback
.. autoclass:: ArtifactSavedEvent
.. autoclass:: ArtifactRemovedEvent
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, Literal, Protocol

import lightning.pytorch as pl
from lightning.pytorch.callbacks.model_checkpoint import ModelCheckpoint
from lightning.pytorch.cli import SaveConfigCallback as _SaveConfigCallback

import lamindb as ln
from lamindb.models.artifact import track_run_input

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
        "mode",
    }
)
_ARTIFACT_AUTO_FEATURES: Final = frozenset(
    {
        "is_best_model",
        "is_last_model",
        "score",
        "model_rank",
        "save_weights_only",
        "monitor",
        "mode",
    }
)
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

    `artifact` is intentionally typed generically so downstream integrations can
    expose their own persisted object while still using the common lifecycle API.
    `storage_uri` is the stable hand-off value for registries such as ClearML.
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

    Artifact-level features:

    - `is_best_model` (bool): Whether this checkpoint is the best model.
    - `is_last_model` (bool): Whether this checkpoint is the most recently saved model.
    - `score` (float): The monitored metric score.
    - `model_rank` (int): Rank among all checkpoints (0 = best).
    - `save_weights_only` (bool): Whether this checkpoint only stores model weights.
    - `monitor` (str): Metric name this checkpoint uses for comparison.
    - `mode` (str): Optimization mode (`min` or `max`) used for checkpoint ranking.

    Run-level features:

    - `logger_name` (str): Name from the first Lightning logger.
    - `logger_version` (str): Version from the first Lightning logger.
    - `max_epochs` (int): Maximum number of epochs.
    - `max_steps` (int): Maximum number of training steps.
    - `precision` (str): Training precision (e.g., "32", "16-mixed", "bf16").
    - `accumulate_grad_batches` (int): Number of batches to accumulate gradients over.
    - `gradient_clip_val` (float): Gradient clipping value.
    - `monitor` (str): Metric name being monitored.
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
    ln.Feature(name="is_last_model", dtype=bool, type=lightning_feature_type).save()
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


class FeatureAnnotator:
    """Manages Lightning feature discovery, collection, and annotation.

    This helper encapsulates all feature-related state and logic used by
    :class:`Checkpoint`.  It handles:

    - Validation of user-specified features at setup time
    - Discovery of auto-features created by :func:`save_lightning_features`
    - Collection of run-level and checkpoint-level feature values
    - Best-model flag management and model rank updates

    The annotator is decoupled from `ModelCheckpoint` state — checkpoint-specific
    values (`best_model_path`, `current_score`, `mode`, etc.) are passed as
    explicit arguments to collection methods.
    """

    def __init__(
        self,
        features: dict[Literal["run", "artifact"], dict[str, Any]] | None = None,
    ) -> None:
        user_features = features or {}
        if invalid_keys := set(user_features) - {"run", "artifact"}:  # type: ignore
            raise ValueError(
                f"Invalid feature keys: {invalid_keys}. Use 'run' and/or 'artifact'."
            )
        self._run_features: dict[str, Any] = user_features.get("run", {})
        self._artifact_features: dict[str, Any] = user_features.get("artifact", {})
        self._auto_features: dict[str, ln.Feature] = {}
        self._hparam_features_available: set[str] = set()
        self._run_features_saved = False

    def setup(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        """Validate user features and discover auto-features.

        Must be called during `Checkpoint.setup()` while `trainer.is_global_zero`
        is `True`.
        """
        self._validate_user_features()
        self._attach_user_run_features()
        self._discover_auto_features()
        self._discover_hparam_features(trainer, pl_module)

    def _attach_user_run_features(self) -> None:
        """Attach user-specified run features to the active LaminDB run."""
        if ln.context.run and self._run_features:
            ln.context.run.features.add_values(self._run_features)

    def _validate_user_features(self) -> None:
        """Ensure all user-specified feature names exist in the database."""
        all_feature_names = set(self._run_features) | set(self._artifact_features)
        if not all_feature_names:
            return
        existing = set(
            ln.Feature.filter(name__in=all_feature_names).values_list("name", flat=True)
        )
        missing = [n for n in all_feature_names if n not in existing]
        if missing:
            s = "s" if len(missing) > 1 else ""
            raise ValueError(
                f"Feature{s} {', '.join(missing)} missing. "
                f"Create {'them' if len(missing) > 1 else 'it'} first."
            )

    def _discover_auto_features(self) -> None:
        """Load auto-features scoped to the `lamindb.lightning` feature type."""
        lightning_feature_type = ln.Feature.filter(
            name="lamindb.lightning", is_type=True
        ).one_or_none()
        self._auto_features.clear()
        if lightning_feature_type is not None:
            self._auto_features = {
                f.name: f
                for f in ln.Feature.filter(
                    name__in=_SUPPORTED_AUTO_FEATURES,
                    type=lightning_feature_type,
                )
            }

    def _discover_hparam_features(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule
    ) -> None:
        """Find which hyperparameter names have matching Features in the DB."""
        hparam_names = self._collect_hparam_names(pl_module, trainer.datamodule)
        self._hparam_features_available = (
            set(ln.Feature.filter(name__in=hparam_names).values_list("name", flat=True))
            if hparam_names
            else set()
        )

    @staticmethod
    def _collect_hparam_names(*sources: Any) -> set[str]:
        """Gather hyperparameter names from one or more sources."""
        names: set[str] = set()
        for source in sources:
            if source is not None and hasattr(source, "hparams") and source.hparams:
                names.update(source.hparams.keys())
        return names

    def get(self, name: str) -> ln.Feature | None:
        """Return the typed auto-feature for *name*, or `None`."""
        return self._auto_features.get(name)

    def _set(self, target: dict[str | ln.Feature, Any], name: str, value: Any) -> None:
        """Add *value* to *target* if the auto-feature *name* is tracked and *value* is not `None`."""
        if (feature := self.get(name)) and value is not None:
            target[feature] = value

    def save_run_features(
        self,
        trainer: pl.Trainer,
        monitor: str | None,
        mode: str,
    ) -> None:
        """Collect and attach run-level features once per run.

        Idempotent — subsequent calls are no-ops.
        """
        if not ln.context.run or self._run_features_saved:
            return

        run_features = self._collect_run_features(trainer, monitor, mode)
        if run_features:
            ln.context.run.features.add_values(run_features)
        self._run_features_saved = True

    def _collect_run_features(
        self,
        trainer: pl.Trainer,
        monitor: str | None,
        mode: str,
    ) -> dict[str | ln.Feature, Any]:
        """Build the dict of run-level feature values (pure, no DB writes)."""
        run_features: dict[str | ln.Feature, Any] = {}

        if trainer.loggers:
            self._set(run_features, "logger_name", trainer.loggers[0].name)
            version = trainer.loggers[0].version
            self._set(
                run_features,
                "logger_version",
                version if isinstance(version, str) else f"version_{version}",
            )

        # Trainer config values
        self._add_trainer_config_features(run_features, trainer, monitor, mode)

        # Hyperparameters
        self._add_hparam_features(
            run_features, trainer.lightning_module, trainer.datamodule
        )

        return run_features

    def _add_trainer_config_features(
        self,
        target: dict[str | ln.Feature, Any],
        trainer: pl.Trainer,
        monitor: str | None,
        mode: str,
    ) -> None:
        """Append trainer configuration values to *target*."""
        self._set(target, "max_epochs", trainer.max_epochs)
        self._set(target, "max_steps", trainer.max_steps)
        self._set(target, "precision", str(trainer.precision))
        self._set(target, "accumulate_grad_batches", trainer.accumulate_grad_batches)
        self._set(target, "gradient_clip_val", trainer.gradient_clip_val)
        self._set(target, "monitor", monitor)
        self._set(target, "mode", mode)

    def _add_hparam_features(
        self,
        target: dict[str | ln.Feature, Any],
        *sources: Any,
    ) -> None:
        """Append hyperparameter values from one or more sources to *target*."""
        for source in sources:
            if source is None:
                continue
            if hasattr(source, "hparams") and source.hparams:
                for name, value in source.hparams.items():
                    if name in self._hparam_features_available:
                        target[name] = value

    def collect_checkpoint_features(
        self,
        trainer: pl.Trainer,
        is_best: bool,
        current_score: Any | None,
        save_weights_only: bool,
        monitor: str | None,
        mode: str,
    ) -> dict[str | ln.Feature, Any]:
        """Collect feature values for a checkpoint artifact.

        All `ModelCheckpoint` state is passed as explicit arguments so the
        annotator stays decoupled from the callback class hierarchy.

        Does **not** mutate existing artifacts — call
        :meth:`clear_best_model_flags` or :meth:`clear_last_model_flags`
        separately when needed.
        """
        feature_values: dict[str | ln.Feature, Any] = {}

        self._set(feature_values, "is_best_model", is_best)
        self._set(feature_values, "is_last_model", True)

        if current_score is not None:
            score = current_score
            if hasattr(score, "item"):
                score = score.item()
            self._set(feature_values, "score", float(score))
        self._set(feature_values, "save_weights_only", save_weights_only)
        self._set(feature_values, "monitor", monitor)
        self._set(feature_values, "mode", mode)

        # User-specified artifact features
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

    def clear_best_model_flags(self, checkpoint_key_prefix: str) -> None:
        """Set `is_best_model=False` on previous best checkpoints."""
        self._clear_flagged_model_feature("is_best_model", checkpoint_key_prefix)

    def clear_last_model_flags(self, checkpoint_key_prefix: str) -> None:
        """Set `is_last_model=False` on previous latest checkpoints."""
        self._clear_flagged_model_feature("is_last_model", checkpoint_key_prefix)

    def _clear_flagged_model_feature(
        self,
        feature_name: Literal["is_best_model", "is_last_model"],
        checkpoint_key_prefix: str,
    ) -> None:
        """Set a boolean model flag to `False` on previously flagged checkpoints."""
        feature = self.get(feature_name)
        if feature is None:
            return
        feature_rows = self._get_artifact_feature_rows(
            {feature_name}, checkpoint_key_prefix
        )
        artifact_ids = [
            artifact_id
            for artifact_id, values in feature_rows.items()
            if values.get(feature_name) is True
        ]
        if not artifact_ids:
            return
        artifacts_by_id = {a.id: a for a in ln.Artifact.filter(id__in=artifact_ids)}
        for artifact_id in artifact_ids:
            if artifact_id not in artifacts_by_id:
                continue
            artifact = artifacts_by_id[artifact_id]
            artifact.features.remove_values(feature, value=True)
            artifact.features.add_values({feature: False})

    def update_model_ranks(self, checkpoint_key_prefix: str, mode: str) -> None:
        """Re-rank all checkpoint artifacts under *checkpoint_key_prefix*."""
        model_rank_feature = self.get("model_rank")
        if model_rank_feature is None:
            return
        feature_rows = self._get_artifact_feature_rows(
            {"score", "model_rank"}, checkpoint_key_prefix
        )
        scored = []
        for artifact_id, values in feature_rows.items():
            if "score" in values:
                scored.append((values["score"], values.get("model_rank"), artifact_id))
        scored.sort(key=lambda x: x[0], reverse=(mode == "max"))

        artifact_ids = [artifact_id for _, _, artifact_id in scored]
        artifacts_by_id = {a.id: a for a in ln.Artifact.filter(id__in=artifact_ids)}
        for rank, (_, old_rank, artifact_id) in enumerate(scored):
            if artifact_id not in artifacts_by_id:
                continue
            af = artifacts_by_id[artifact_id]
            if old_rank is not None:
                af.features.remove_values(model_rank_feature, value=old_rank)
            af.features.add_values({model_rank_feature: rank})

    def _get_artifact_feature_rows(
        self,
        feature_names: set[str],
        checkpoint_key_prefix: str,
    ) -> dict[int, dict[str, Any]]:
        """Query feature values for checkpoint artifacts under *checkpoint_key_prefix*.

        Returns a dict keyed by artifact ID, where each value is a dict mapping
        feature name to its stored value.  Example::

            {
                42: {"score": 0.95, "is_best_model": True},
                71: {"score": 0.87, "is_best_model": False, "model_rank": 1},
            }
        """
        feature_ids = [
            feature.id for name in feature_names if (feature := self.get(name))
        ]
        key_startswith = checkpoint_key_prefix + "/"
        if feature_ids:
            rows = ln.models.ArtifactJsonValue.filter(
                artifact__key__startswith=key_startswith,
                jsonvalue__feature_id__in=feature_ids,
            ).values_list("artifact_id", "jsonvalue__feature__name", "jsonvalue__value")
        else:
            rows = ln.models.ArtifactJsonValue.filter(
                artifact__key__startswith=key_startswith,
                jsonvalue__feature__name__in=feature_names,
            ).values_list("artifact_id", "jsonvalue__feature__name", "jsonvalue__value")
        result: dict[int, dict[str, Any]] = {}
        for artifact_id, feature_name, value in rows:
            if artifact_id not in result:
                result[artifact_id] = {}
            result[artifact_id][feature_name] = value
        return result


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

    def _notify_artifact_saved(
        self,
        trainer: pl.Trainer,
        *,
        kind: ArtifactKind,
        key: str,
        artifact: Any,
        local_path: Path | str,
    ) -> ArtifactSavedEvent:
        event = ArtifactSavedEvent(
            kind=kind,
            key=key,
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
    explicit `dirpath` or the trainer's logger configuration.

    When `dirpath` is omitted (recommended), Lightning decides where to store
    checkpoints locally (typically `lightning_logs/version_N/checkpoints/`)
    and the artifact key is derived from the logger's `save_dir`, `name`,
    and `version`.  When `dirpath` is provided, it is used directly as the
    key prefix.

    All artifacts are scoped under a single **base prefix**.  Checkpoints
    (and `hparams.yaml`) live under `{base}/checkpoints/`; other artifacts
    (e.g. `config.yaml`) live directly under `{base}/`.

    Base prefix derivation (highest priority first):

    1. `dirpath` provided → `{dirpath}` (logger is ignored for key purposes)
    2. `dirpath` omitted, logger present → `{save_dir_basename}/{name}/{version}`
    3. `dirpath` omitted, no logger → empty

    When `run_uid_is_version` is `True` (the default) and a Lamin run context
    is active, the run UID is incorporated into the base prefix:

    - Case 1/3: the run UID is appended as an extra path segment
      (e.g. `my/dir/{run_uid}`, or just `{run_uid}`).
    - Case 2: the logger's auto-incremented `version` is *replaced* by the
      run UID (`{save_dir_basename}/{name}/{run_uid}`).

    Resulting key layout (with run UID active)::

        {base}/checkpoints/epoch=0-step=100.ckpt
        {base}/checkpoints/hparams.yaml
        {base}/config.yaml

    If available in the database through `save_lightning_features()`, the following `lamindb.lightning` features are automatically tracked:

    - Artifact-level: `is_best_model`, `is_last_model`, `score`, `model_rank`, `save_weights_only`, `monitor`, `mode`
    - Run-level: `logger_name`, `logger_version`, `max_epochs`, `max_steps`, `precision`, `accumulate_grad_batches`, `gradient_clip_val`, `monitor`, `mode`

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
        run_uid_is_version: When `True` (default) and a Lamin run context is
            active, incorporate the run UID into the base prefix. For the
            logger case the logger's auto-incremented version is replaced;
            for the dirpath and no-logger cases the run UID is appended as
            an extra path segment. Prevents cross-run key collisions.
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

        Explicit `dirpath` for full control over the artifact key prefix::

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

        For more, see the guide: :doc:`lightning`.
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
        run_uid_is_version: bool = True,
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
        self._feature_annotator = FeatureAnnotator(features)
        self._hparams_yaml_saved = False
        self._run_uid_is_version = run_uid_is_version
        self._trainer: pl.Trainer | None = None
        self._artifact_publisher: ArtifactPublisher = LaminArtifactPublisher()

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        """Validate user features and detect available auto-features."""
        super().setup(trainer, pl_module, stage)
        self._trainer = trainer

        if self.save_last:
            warnings.warn(
                "save_last is not necessary with Lamin. Checkpoint metadata"
                " (is_best_model, is_last_model, model_rank, score) makes the latest checkpoint"
                " queryable without encoding this in the filename. Consider"
                " disabling save_last to avoid redundant checkpoint copies.",
                UserWarning,
                stacklevel=2,
            )

        if trainer.is_global_zero:
            self._feature_annotator.setup(trainer, pl_module)

    def _base_prefix(self, trainer: pl.Trainer) -> str:
        """Compute the base artifact key prefix.

        The base prefix is the root namespace for all artifacts produced by
        this callback.  Checkpoints live under `{base}/checkpoints/` and
        other files (config, hparams) directly under `{base}/`.

        Priority: explicit `dirpath` > logger > run UID > empty.
        """
        run_uid = self._active_run_uid()
        if self._original_dirpath is not None:
            prefix = str(self._original_dirpath).rstrip("/")
            return f"{prefix}/{run_uid}" if run_uid else prefix
        if len(trainer.loggers) > 0:
            return self._logger_prefix(trainer, run_uid)
        return run_uid or ""

    def _active_run_uid(self) -> str | None:
        """Return the Lamin run UID when run-UID scoping is active."""
        if self._run_uid_is_version and ln.context.run is not None:
            return ln.context.run.uid
        return None

    def _logger_prefix(self, trainer: pl.Trainer, run_uid: str | None) -> str:
        """Derive a key prefix from the trainer's first logger."""
        assert trainer.loggers, "_logger_prefix requires at least one logger"
        logger = trainer.loggers[0]
        save_dir = logger.save_dir or trainer.default_root_dir
        name = str(logger.name).rstrip("/")
        if run_uid:
            version = run_uid
        else:
            version = logger.version
            version = version if isinstance(version, str) else f"version_{version}"
        return f"{Path(save_dir).name}/{name}/{version.rstrip('/')}"

    @property
    def base_prefix(self) -> str:
        """The base artifact key prefix for all artifacts from this callback.

        Checkpoints live under `{base_prefix}/checkpoints/` and configs
        directly under `{base_prefix}/`.

        Available after `setup()` has been called.
        """
        assert self._trainer is not None, "base_prefix is only available after setup()"
        return self._base_prefix(self._trainer)

    @property
    def checkpoint_key_prefix(self) -> str:
        """The artifact key prefix used for checkpoint artifacts.

        Available after `setup()` has been called, for example once
        `trainer.fit()` has started.
        """
        base = self.base_prefix
        return f"{base}/checkpoints" if base else "checkpoints"

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
        base = self._base_prefix(trainer)
        if kind in {"checkpoint", "hparams"}:
            prefix = f"{base}/checkpoints" if base else "checkpoints"
        else:
            prefix = base
        if prefix:
            return f"{prefix}/{Path(filepath).name}"
        return Path(filepath).name

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
        self._feature_annotator.clear_last_model_flags(self.checkpoint_key_prefix)

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
        key = self.resolve_artifact_key(
            trainer=trainer, filepath=filepath, kind="checkpoint"
        )
        existing_artifact = ln.Artifact.filter(key=key).one_or_none()
        if existing_artifact is not None:
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
            key=key,
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
        key = self.resolve_artifact_key(
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
            key=key,
            artifact=artifact,
            local_path=config_path,
        )
        return artifact

    def save_hparams_artifact(
        self, trainer: pl.Trainer, hparams_path: Path | str
    ) -> ln.Artifact | None:
        """Save Lightning's auto-generated hparams file and emit the event.

        Returns `None` if Lightning did not generate `hparams.yaml` for the
        current run.
        """
        if not Path(hparams_path).exists():
            return None

        key = self.resolve_artifact_key(
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
            key=key,
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

    def _save_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Save checkpoint to the instance."""
        super()._save_checkpoint(trainer, filepath)

        if not trainer.is_global_zero:
            return
        self._save_hparams_yaml(trainer)

        self._feature_annotator.save_run_features(
            trainer, monitor=self.monitor, mode=self.mode
        )
        self._feature_annotator.clear_last_model_flags(self.checkpoint_key_prefix)
        is_best = self.best_model_path == str(filepath)
        feature_values = self._feature_annotator.collect_checkpoint_features(
            trainer,
            is_best=is_best,
            current_score=self.current_score,
            save_weights_only=self.save_weights_only,
            monitor=self.monitor,
            mode=self.mode,
        )

        if is_best:
            self._feature_annotator.clear_best_model_flags(self.checkpoint_key_prefix)

        self.save_checkpoint_artifact(trainer, filepath, feature_values=feature_values)

        self._feature_annotator.update_model_ranks(
            self.checkpoint_key_prefix, mode=self.mode
        )

    def _remove_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Remove the local checkpoint file and emit a removal event."""
        artifact: ln.Artifact | None = None
        key = self.resolve_artifact_key(
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
            if artifact is not None:
                artifact.delete(permanent=True, storage=True)


class SaveConfigCallback(_SaveConfigCallback):
    """SaveConfigCallback that also saves config to the instance.

    Use with LightningCLI to save the resolved configuration file alongside checkpoints.

    The local config file is saved under `{save_dir}/{name}/{version}/`
    derived from the first logger, avoiding Lightning's `trainer.log_dir`
    which hardcodes an `isinstance` check for `TensorBoardLogger` /
    `CSVLogger` and silently changes the directory for other loggers.

    This callback looks for any :class:`ArtifactPublishingModelCheckpoint`, not just
    Lamin's concrete :class:`Checkpoint`. That keeps the config-save path aligned
    with custom subclasses built on the generic artifact-publishing base.

    Config artifacts are stored directly under the **base prefix** of the
    active :class:`Checkpoint` callback.  The base prefix follows the same
    derivation rules as for checkpoints (dirpath > logger > empty), so
    configs are always co-located with their checkpoints:

    - `Checkpoint.dirpath` set → `{dirpath}/config.yaml`
      (`{dirpath}/{run_uid}/config.yaml` with run-UID scoping)
    - Logger present, no `dirpath` → `{save_dir_basename}/{name}/{version}/config.yaml`
    - Neither → `config.yaml` (or `{run_uid}/config.yaml` with run-UID scoping)

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

        if self.save_to_log_dir:
            config_path = self._config_path(trainer)

            if not self.overwrite:
                file_exists = config_path.exists() if trainer.is_global_zero else False
                file_exists = trainer.strategy.broadcast(file_exists)
                if file_exists:
                    raise RuntimeError(f"Config file already exists: {config_path}")

            if trainer.is_global_zero:
                config_path.parent.mkdir(exist_ok=True, parents=True)
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

    def _config_path(self, trainer: pl.Trainer) -> Path:
        """Derive the local config file path from the first logger.

        We intentionally avoid `trainer.log_dir` because Lightning hardcodes
        an `isinstance` check against `TensorBoardLogger` and `CSVLogger`
        there.  For those two loggers it uses `logger.log_dir` (which appends
        name/version), while for every other logger it falls back to
        `logger.save_dir` (no name/version).  This means the config file
        location silently changes depending on which logger happens to be first
        — making it unpredictable for third-party loggers.

        This method always uses `logger.save_dir` + `name` + `version`,
        giving a consistent directory layout regardless of logger type.
        """
        if len(trainer.loggers) > 0:
            first = trainer.loggers[0]
            save_dir = (
                first.save_dir
                if first.save_dir is not None
                else trainer.default_root_dir
            )
            name = first.name
            version = first.version
            version = version if isinstance(version, str) else f"version_{version}"
            return Path(save_dir) / str(name) / version / self.config_filename
        return Path(trainer.default_root_dir) / self.config_filename

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
