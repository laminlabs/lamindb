from pathlib import Path
from typing import Optional

import lamindb as ln
import lightning.pytorch as pl
import torch
from lightning.fabric.utilities.types import _PATH
from lightning.pytorch.callbacks.model_checkpoint import ModelCheckpoint

from lamindb.integrations.lightning.util import track_if_not_tracked


class LaminCheckpoint(ModelCheckpoint):
    r"""
    Pytorch Lightning checkpointing callback that saves transfers the checkpoint to LaminDB and
    adds some information, like the current score, as lamin "features".
    Because Lamin is versioning files, it is not necessary to save the checkpoints under
    different names while the training progresses.

    When you set **save_top_k** to a value other than ``-1``, the best k checkpoints are
    kept in the lamin version history and the rest are pruned.

    The latest checkpoint is also kept, independent of whether it is among the best or not.

    Save the model periodically by monitoring a quantity. Every metric logged with
    :meth:`~lightning.pytorch.core.LightningModule.log` or :meth:`~lightning.pytorch.core.LightningModule.log_dict` is
    a candidate for the monitor key. For more information, see :ref:`checkpointing`.

    After training finishes, use :attr:`best_model_path` to retrieve the path to the
    best checkpoint file and :attr:`best_model_score` to retrieve its score.

    Args:
        lamin_instance: LaminDB instance to connect to (example: "pfizer/ib-mlr").
        lamin_root_dir: Directory to save the model checkpoints in LaminDB (example: "my-model-logs").
        lamin_file_name: Name of the file to save the model checkpoint as in LaminDB (example: "model.ckpt").
        min_epochs: checkpoints for runs with less epochs than `min_epochs` will be removed from lamin.
        dirpath: directory to locally save the model file.

            Example::

                # custom path
                # saves a file like: my/path/epoch=0-step=10.ckpt
                >>> checkpoint_callback = ModelCheckpoint(
                ...                         lamin_instance="pfizer/ib-mlr",  # Example LaminDB instance
                ...                         lamin_root_dir="example_experiment",
                ...                         lamin_file_name = "model.ckpt"
                ...                         dirpath='my/path/'
                ...                     )
            By default, dirpath is ``None`` and will be set at runtime to the location
            specified by :class:`~lightning.pytorch.trainer.trainer.Trainer`'s
            :paramref:`~lightning.pytorch.trainer.trainer.Trainer.default_root_dir` argument,
            and if the Trainer uses a logger, the path will also contain logger name and version.

        filename: checkpoint filename. Can contain named formatting options to be auto-filled.

            Example::

                # save any arbitrary metrics like `val_loss`, etc. in name
                # saves a file like: my/path/epoch=2-val_loss=0.02-other_metric=0.03.ckpt
                >>> checkpoint_callback = ModelCheckpoint(
                ...                         lamin_instance="pfizer/ib-mlr",  # Example LaminDB instance
                ...                         lamin_root_dir="example_experiment",
                ...                         lamin_file_name = "model.ckpt"
                ...                         filename='{epoch}-{val_loss:.2f}-{other_metric:.2f}'
                ...                     )
            By default, filename is ``None`` and will be set to ``'{epoch}-{step}'``, where "epoch" and "step" match
            the number of finished epoch and optimizer steps respectively.
        monitor: quantity to monitor. By default it is ``None`` which saves a checkpoint only for the last epoch.
        verbose: verbosity mode. Default: ``False``.
        save_last: When ``True``, saves a `last.ckpt` copy whenever a checkpoint file gets saved. Can be set to
            ``'link'`` on a local filesystem to create a symbolic link. This allows accessing the latest checkpoint
            in a deterministic manner. Default: ``None``.
        save_top_k: if ``save_top_k == k``,
            the best k models according to the quantity monitored will be saved.
            If ``save_top_k == 0``, no models are saved.
            If ``save_top_k == -1``, all models are saved.
            Please note that the monitors are checked every ``every_n_epochs`` epochs.
            If ``save_top_k >= 2`` and the callback is called multiple times inside an epoch, and the filename remains
            unchanged, the name of the saved file will be appended with a version count starting with ``v1`` to avoid
            collisions unless ``enable_version_counter`` is set to False. The version counter is unrelated to the top-k
            ranking of the checkpoint, and we recommend formatting the filename to include the monitored metric to avoid
            collisions.
        mode: one of {min, max}.
            If ``save_top_k != 0``, the decision to overwrite the current save file is made
            based on either the maximization or the minimization of the monitored quantity.
            For ``'val_acc'``, this should be ``'max'``, for ``'val_loss'`` this should be ``'min'``, etc.
        auto_insert_metric_name: When ``True``, the checkpoints filenames will contain the metric name.
            For example, ``filename='checkpoint_{epoch:02d}-{acc:02.0f}`` with epoch ``1`` and acc ``1.12`` will resolve
            to ``checkpoint_epoch=01-acc=01.ckpt``. Is useful to set it to ``False`` when metric names contain ``/``
            as this will result in extra folders.
            For example, ``filename='epoch={epoch}-step={step}-val_acc={val/acc:.2f}', auto_insert_metric_name=False``
        save_weights_only: if ``True``, then only the model's weights will be
            saved. Otherwise, the optimizer states, lr-scheduler states, etc are added in the checkpoint too.
        every_n_train_steps: Number of training steps between checkpoints.
            If ``every_n_train_steps == None or every_n_train_steps == 0``, we skip saving during training.
            To disable, set ``every_n_train_steps = 0``. This value must be ``None`` or non-negative.
            This must be mutually exclusive with ``train_time_interval`` and ``every_n_epochs``.
        train_time_interval: Checkpoints are monitored at the specified time interval.
            For all practical purposes, this cannot be smaller than the amount
            of time it takes to process a single training batch. This is not
            guaranteed to execute at the exact time specified, but should be close.
            This must be mutually exclusive with ``every_n_train_steps`` and ``every_n_epochs``.
        every_n_epochs: Number of epochs between checkpoints.
            This value must be ``None`` or non-negative.
            To disable saving top-k checkpoints, set ``every_n_epochs = 0``.
            This argument does not impact the saving of ``save_last=True`` checkpoints.
            If all of ``every_n_epochs``, ``every_n_train_steps`` and
            ``train_time_interval`` are ``None``, we save a checkpoint at the end of every epoch
            (equivalent to ``every_n_epochs = 1``).
            If ``every_n_epochs == None`` and either ``every_n_train_steps != None`` or ``train_time_interval != None``,
            saving at the end of each epoch is disabled
            (equivalent to ``every_n_epochs = 0``).
            This must be mutually exclusive with ``every_n_train_steps`` and ``train_time_interval``.
            Setting both ``ModelCheckpoint(..., every_n_epochs=V, save_on_train_epoch_end=False)`` and
            ``Trainer(max_epochs=N, check_val_every_n_epoch=M)``
            will only save checkpoints at epochs 0 < E <= N
            where both values for ``every_n_epochs`` and ``check_val_every_n_epoch`` evenly divide E.
        save_on_train_epoch_end: Whether to run checkpointing at the end of the training epoch.
            If this is ``False``, then the check runs at the end of the validation.
        enable_version_counter: Whether to append a version to the existing file name.
            If this is ``False``, then the checkpoint files will be overwritten.

    Note:
        For extra customization, ModelCheckpoint includes the following attributes:

        - ``CHECKPOINT_JOIN_CHAR = "-"``
        - ``CHECKPOINT_EQUALS_CHAR = "="``
        - ``CHECKPOINT_NAME_LAST = "last"``
        - ``FILE_EXTENSION = ".ckpt"``
        - ``STARTING_VERSION = 1``

        For example, you can change the default last checkpoint name by doing
        ``checkpoint_callback.CHECKPOINT_NAME_LAST = "{epoch}-last"``

        If you want to checkpoint every N hours, every M train batches, and/or every K val epochs,
        then you should create multiple ``ModelCheckpoint`` callbacks.

        If the checkpoint's ``dirpath`` changed from what it was before while resuming the training,
        only ``best_model_path`` will be reloaded and a warning will be issued.

    Raises:
        MisconfigurationException:
            If ``save_top_k`` is smaller than ``-1``,
            if ``monitor`` is ``None`` and ``save_top_k`` is none of ``None``, ``-1``, and ``0``, or
            if ``mode`` is none of ``"min"`` or ``"max"``.
        ValueError:
            If ``trainer.save_checkpoint`` is ``None``.

    Example::

        >>> from lightning.pytorch import Trainer
        >>> from lightning.pytorch.callbacks import ModelCheckpoint

        # saves checkpoints to 'my/path/' at every epoch
        >>> checkpoint_callback = ModelCheckpoint(dirpath='my/path/')
        >>> trainer = Trainer(callbacks=[checkpoint_callback])

        # save epoch and val_loss in name
        # saves a file like: my/path/sample-mnist-epoch=02-val_loss=0.32.ckpt
        >>> checkpoint_callback = ModelCheckpoint(
        ...     monitor='val_loss',
        ...     dirpath='my/path/',
        ...     filename='sample-mnist-{epoch:02d}-{val_loss:.2f}'
        ... )

        # save epoch and val_loss in name, but specify the formatting yourself (e.g. to avoid problems with Tensorboard
        # or Neptune, due to the presence of characters like '=' or '/')
        # saves a file like: my/path/sample-mnist-epoch02-val_loss0.32.ckpt
        >>> checkpoint_callback = ModelCheckpoint(
        ...     monitor='val/loss',
        ...     dirpath='my/path/',
        ...     filename='sample-mnist-epoch{epoch:02d}-val_loss{val/loss:.2f}',
        ...     auto_insert_metric_name=False
        ... )

        # retrieve the best checkpoint after training
        checkpoint_callback = ModelCheckpoint(dirpath='my/path/')
        trainer = Trainer(callbacks=[checkpoint_callback])
        model = ...
        trainer.fit(model)
        checkpoint_callback.best_model_path

    .. tip:: Saving and restoring multiple checkpoint callbacks at the same time is supported under variation in the
        following arguments:

        *monitor, mode, every_n_train_steps, every_n_epochs, train_time_interval*

        Read more: :ref:`Persisting Callback State <extensions/callbacks_state:save callback state>`

    """

    def __init__(
        self,
        lamin_instance: str,
        lamin_root_dir: Path | str,
        lamin_file_name: str = "model.ckpt",
        min_epochs: int = 1,
        dirpath: Optional[_PATH] = None,
        *args,
        **kwargs,
    ):
        self._dir_path_is_user_defined = dirpath is not None
        super().__init__(dirpath=dirpath, *args, **kwargs)
        self.lamin_instance = lamin_instance
        self.lamin_root_dir = lamin_root_dir
        self.lamin_file_name = lamin_file_name
        self.min_epochs = min_epochs

        self._logger_labels = {}
        self._created_features = set()

    @property
    def _lamin_checkpoint_dir(self):
        """
        Returns the directory in LaminDB where the checkpoints are stored. What this
        should be depends on whether the user activity defined a 'dirpath` or not.

        """
        dirpath = Path(self.dirpath)

        if self._dir_path_is_user_defined:
            # dirpath is user defined, deciding how to "mirror"
            # on lamin is more ambiguous.
            if dirpath.is_absolute():
                if dirpath.is_relative_to(Path.cwd()):
                    # if the dirpath is absolute but relative to the current
                    # working directory, we interpret it as the path relative
                    # to the current working directory.
                    return Path(self.lamin_root_dir) / dirpath.relative_to(Path.cwd())
                else:
                    # if the dirpath is absolute and not relative to the current
                    #  working directory, we don't know what to do, so we just use
                    #  the lamin_root_dir.
                    return Path(self.lamin_root_dir)
            else:
                # if the dirpath is a relative path, we interpret it as
                # relative to the lamin_root_dir
                return Path(self.lamin_root_dir) / dirpath

        elif len(self._logger_labels) > 0:
            # checkpoint dir is inherited from the first logger's
            # .name and .version
            return Path(self.lamin_root_dir) / dirpath.relative_to(
                dirpath.parent.parent.parent
            )
        else:
            # no logger, checkpoints end up in trainer.default_root_dir/checkpoints
            # which we mirror as "lamin_root_dir/checkpoints"
            return Path(self.lamin_root_dir) / "checkpoints"

    @property
    def _lamin_path(self):
        return Path(self._lamin_checkpoint_dir) / self.lamin_file_name

    def setup(
        self, trainer: pl.Trainer, pl_module: pl.LightningModule, stage: str
    ) -> None:
        # Connect to LaminDB

        super().setup(trainer=trainer, pl_module=pl_module, stage=stage)

        self.setup_logger_labels(trainer=trainer)
        if trainer.is_global_zero:
            track_if_not_tracked(self.lamin_instance)
            # if not self._run:
            #     ln.connect(self.lamin_instance)
            #     transform = ln.Transform("model_training")
            #     transform.save()
            #     self._run = ln.Run(transform=transform)
            #     self._run.save()
            self._setup_lamin_features()
        print(f"Transferring checkpoints to LaminDB at {self._lamin_path}")

    def setup_logger_labels(self, trainer: pl.Trainer) -> None:
        for logger in trainer.loggers:
            class_name = logger.__class__.__name__
            self._logger_labels[f"{class_name}.name"] = logger.name
            # version = logger.version if isinstance(logger.version, str) else f"version_{logger.version}"
            self._logger_labels[f"{class_name}.version"] = logger.version

    def _setup_lamin_features(self):
        ln.Feature(name="is_best_model", dtype=bool, default_value=False).save()
        ln.Feature(name="score", dtype=float).save()
        ln.Feature(name="model_rank", dtype=int).save()

        for key, value in self._logger_labels.items():
            ln.Feature(name=key, dtype=type(value)).save()

    def _assure_feature(self, name, dtype) -> None:
        if name in self._created_features:
            return
        ln.Feature(name=name, dtype=dtype).save()
        self._created_features.add(name)

    def _save_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        # Save the checkpoint using the parent class
        super()._save_checkpoint(trainer, filepath)
        self._lamin_upload_last(trainer, filepath)

    def _lamin_upload_last(self, trainer: pl.Trainer, filepath: str) -> None:
        if not trainer.is_global_zero:
            return
        artifact = ln.Artifact(
            filepath, key=str(self._lamin_path), description="Latest Model checkpoint"
        )
        artifact.save()

        artifact.features.add_values(self._logger_labels)

        if self.best_model_path == filepath:
            self._erase_best_model_label()
            artifact.features.add_values({"is_best_model": True})

        monitor_candidates = {}
        for key, value in self._monitor_candidates(trainer).items():
            if torch.is_tensor(value):
                value = value.item() if value.numel() == 1 else value.tolist()
            dtype = type(value)
            if dtype is str:
                dtype = ln.ULabel
            monitor_candidates[key] = value
            self._assure_feature(name=key, dtype=dtype)
        artifact.features.add_values(monitor_candidates)

        artifact.features.add_values({"score": self.current_score.item()})

        self._prune_checkpoint_history()

    def _erase_best_model_label(self):
        """
        Goes into the existing checkpoints in the lamin DB and
        sets the `is_best_model` feature to False.
        """

        for artifact in self.get_best_model_artifacts():
            artifact.features.remove_values("is_best_model")

    def _prune_checkpoint_history(self):
        prune_artifacts = self.save_top_k != -1

        artifacts = ln.Artifact.filter(
            key=str(self._lamin_path)
        )  # .filter(is_latest=False)
        sorted_atrifacts = sorted(
            [
                (artifact.features.get_values()["score"], artifact)
                for artifact in artifacts
            ],
            reverse=self.mode == "max",
        )

        best_scoring = sorted_atrifacts
        if prune_artifacts:
            best_scoring = sorted_atrifacts[: self.save_top_k]

        for rank, (_, artifact) in enumerate(best_scoring):
            artifact.features.remove_values("model_rank")
            artifact.features.add_values({"model_rank": rank})

        if not prune_artifacts:
            return

        for _, artifact in sorted_atrifacts[self.save_top_k :]:
            if artifact.is_latest:
                continue
            artifact.delete(permanent=True, storage=True)

    def get_best_model_artifacts(self):
        """This should in principle just return a set with one element."""

        return ln.Artifact.filter(is_best_model=True).filter(key=str(self._lamin_path))

    def on_exception(
        self,
        trainer: "pl.Trainer",
        pl_module: "pl.LightningModule",
        exception: BaseException,
    ) -> None:
        self.delete_short_run(trainer=trainer)

    def on_fit_end(
        self, trainer: "pl.Trainer", pl_module: "pl.LightningModule"
    ) -> None:
        self.delete_short_run(trainer=trainer)

    def delete_short_run(self, trainer: pl.Trainer) -> None:
        """
        Deletes the run if it is shorter than self.min_epochs.
        This is useful to avoid cluttering the LaminDB with runs
        that were not actually trained.
        """
        if not trainer.is_global_zero:
            return
        if trainer.current_epoch >= self.min_epochs:
            return
        if not ln.context.run:
            return
        print(f"Deleting short run {ln.context.run.id} from LaminDB")
        artifacts = ln.Artifact.filter(key=str(self._lamin_path))
        for artifact in artifacts:
            print("Deleting Artifact:", artifact.id, "key ", self._lamin_path)
            artifact.delete(permanent=True, storage=True)
        ln.context.run.delete()
