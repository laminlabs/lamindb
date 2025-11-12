import os
from pathlib import Path

import lamindb as ln
from lightning import LightningModule, Trainer
from lightning.fabric.utilities.cloud_io import get_filesystem
from lightning.pytorch.cli import SaveConfigCallback

from lamindb.integrations.lightning.checkpoint import LaminCheckpoint
from lamindb.integrations.lightning.util import track_if_not_tracked


class LaminSaveConfigCallback(SaveConfigCallback):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setup(self, trainer: Trainer, pl_module: LightningModule, stage: str) -> None:
        if self.already_saved:
            return

        if self.save_to_log_dir:
            log_dir = trainer.log_dir  # this broadcasts the directory
            assert log_dir is not None
            config_path = os.path.join(log_dir, self.config_filename)
            fs = get_filesystem(log_dir)

            if not self.overwrite:
                # check if the file exists on rank 0
                file_exists = (
                    fs.isfile(config_path) if trainer.is_global_zero else False
                )
                # broadcast whether to fail to all ranks
                file_exists = trainer.strategy.broadcast(file_exists)
                if file_exists:
                    raise RuntimeError("Config file already exists")
            if trainer.is_global_zero:
                # save only on rank zero to avoid race conditions.
                # the `log_dir` needs to be created as we rely on
                # the logger to do it usually but it hasn't logged
                # anything at this point
                fs.makedirs(log_dir, exist_ok=True)
                self.parser.save(
                    self.config,
                    config_path,
                    skip_none=False,
                    overwrite=self.overwrite,
                    multifile=self.multifile,
                )
                lamin_experiment_dir = self.lamin_experiment_dir(trainer)
                lamin_config_path = lamin_experiment_dir / self.config_filename
                track_if_not_tracked(self.lamin_instance(trainer))
                artifact = ln.Artifact(
                    config_path,
                    key=str(lamin_config_path),
                    description="Config file for lightning CLI",
                )
                artifact.save()

        if trainer.is_global_zero:
            self.save_config(trainer, pl_module, stage)
            self.already_saved = True

        # broadcast so that all ranks are in sync on future calls to .setup()
        self.already_saved = trainer.strategy.broadcast(self.already_saved)

    def get_lamin_checkpoint_callback(self, trainer: Trainer) -> LaminCheckpoint:
        checkpoint_callbacks = [
            cb for cb in trainer.callbacks if isinstance(cb, LaminCheckpoint)
        ]
        if not checkpoint_callbacks:
            raise RuntimeError(
                "LaminSaveConfigCallback requires a LaminCheckpoint callback to be used."
            )

        return checkpoint_callbacks[0]

    def lamin_instance(self, trainer: Trainer) -> str:
        lamin_checkpoint_callback = self.get_lamin_checkpoint_callback(trainer)
        return lamin_checkpoint_callback.lamin_instance

    def lamin_experiment_dir(self, trainer: Trainer) -> Path:
        lamin_checkpoint_callback = self.get_lamin_checkpoint_callback(trainer)

        if lamin_checkpoint_callback._dir_path_is_user_defined:
            dirpath = Path(lamin_checkpoint_callback.dirpath)
            # dirpath is user defined, deciding how to "mirror"
            # on lamin is more ambiguous.
            if dirpath.is_absolute():
                if dirpath.is_relative_to(Path.cwd()):
                    # if the dirpath is absolute but relative to the current
                    # working directory, we interpret it as the path relative
                    # to the current working directory.
                    return Path(
                        lamin_checkpoint_callback.lamin_root_dir
                    ) / dirpath.relative_to(Path.cwd())
                else:
                    # if the dirpath is absolute and not relative to the current
                    #  working directory, we don't know what to do, so we just use
                    #  the lamin_root_dir.
                    return Path(lamin_checkpoint_callback.lamin_root_dir)
            else:
                # if the dirpath is a relative path, we interpret it as
                # relative to the lamin_root_dir
                return Path(lamin_checkpoint_callback.lamin_root_dir) / dirpath

        elif len(trainer.loggers) > 0:
            first_logger = trainer.loggers[0]
            version = first_logger.version
            version = version if isinstance(version, str) else f"version_{version}"
            return (
                Path(lamin_checkpoint_callback.lamin_root_dir)
                / first_logger.name
                / version
            )
        else:
            return Path(lamin_checkpoint_callback.lamin_root_dir)
