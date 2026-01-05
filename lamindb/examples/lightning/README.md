# lamin-lightning

A Python package that integrates [LaminDB](https://lamin.ai) with [PyTorch Lightning](https://lightning.ai/pytorch-lightning) to streamline model checkpointing.

## tl;dr
Take a look at the [example](lamindb/examples/lightning) directory for a full example of how
to use this package with PyTorch Lightning CLI.

In a environment with pytorch, lightning and lamindb installed:

```bash
cd lamindb/examples/lightning
python main.py fit --config lightning_config.yml
```

## Quick Start

```python
from lightning import Trainer
from lamindb.integrations.lightning.checkpoint import LaminCheckpoint

checkpoint_callback = LaminCheckpoint(
    lamin_instance="your/instance",
    lamin_root_dir="checkpoints",
    monitor="val_loss",
    save_top_k=3
)

trainer = Trainer(callbacks=[checkpoint_callback])
```

## Using LaminCheckpoint

The `LaminCheckpoint` extends Lightning's `ModelCheckpoint` to save checkpoints to your LaminDB instance.

### Configuration Options

Configuration options in addition to those provided by `ModelCheckpoint`:

| Option            | Description                                          | Default      |
|-------------------|------------------------------------------------------|--------------|
| `lamin_instance`  | LaminDB instance identifier (e.g., "org/project")    | Required     |
| `lamin_root_dir`  | Root directory for checkpoint storage in LaminDB     | Required     |
| `lamin_file_name` | File name used in Lamin                              | "model.ckpt" |
| `min_epochs`      | Runs with less epochs than min_epochs will be pruned | 0            |

The same lamin_file_name will be used as the key in lamin for all checkpoints of a training run.
Because lamin has versioning of files and metadata labels attached to these individual versions
(like epoch, step and score) we don't have to squeeze this information into the file name like
is done for local checkpoints.

When `save_top_k` is set, the callback will erase old checkpoints from lamin that are not part
of the top k any longer.

### YAML Configuration Example

This is a snippet from a lightning_config.yml (part shown) that can be used by the lightning CLI:

```yaml
trainer:
  callbacks:
    - class_path: lamindb.integrations.lightning.LaminCheckpoint
      init_args:
        lamin_instance: "your/instance"
        lamin_root_dir: "my-checkpoints"
        every_n_epochs: 1
        min_epochs: 2
        save_top_k: 3
        monitor: val_loss
```
together with a script containing something like:

```python
# main.py
from lightning.pytorch.cli import ArgsType, LightningCLI
from lamindb.integrations.lightning.save_config_callback import LaminSaveConfigCallback

def my_cli(args: ArgsType = None):
    LightningCLI(
        model_class=MyLightningModule,
        datamodule_class=MyDataModule,
        save_config_callback=LaminSaveConfigCallback,
        args=args,
    )

if __name__ == "__main__":
    my_cli()

```
Using the `LaminSaveConfigCallback` will automatically save a config.yml to LaminDB
when the training run starts.

For a full example, see the [lightning_config.yml](lightning_config.yml) in the example directory.
This can be used like:

```bash
cd lamindb/examples/lightning
python main.py fit --config lightning_config.yml
```

## Features

- Automatically stores checkpoints in LaminDB
- Logger metadata: for each logger, experiment name and version
  are saved as metadata to the checkpoint in LaminDB
- Supports `save_top_k` to keep only the best checkpoints, also
  in lamin.
- Short runs that don't meet the `min_epochs` requirement
  will be pruned from LaminDB to avoid cluttering
