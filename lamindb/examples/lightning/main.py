from typing import Dict, Iterator

import lightning.pytorch as pl
import torch
from lightning.pytorch import LightningDataModule
from lightning.pytorch.cli import ArgsType, LightningCLI
from torch.utils.data import DataLoader, IterableDataset

from lamin_lightning import LaminSaveConfigCallback


class ExampleDataset(IterableDataset):
    """An example pytorch Dataset object returning
       some spoof inputs and targets based on the function
       ax**2 + b. The task of the neural network will be to
       predict whether a is positive or negative.
       The task of the neural network will be to
    Args:
        balance: change for training example to
        be of the positive class. Defaults to 0.5 (balanced).

    """

    def __init__(self, balance: float = 0.5, batch_size: int = 64, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.balance = balance
        self.batch_size = batch_size

    def __len__(self):
        return 100

    def __iter__(self) -> Iterator[Dict[str, torch.Tensor]]:
        self._n = 0
        return self

    def __next__(self):
        if self._n >= len(self):
            raise StopIteration

        is_positive = torch.bernoulli(torch.tensor([self.balance] * self.batch_size))
        sign = is_positive * 2 - 1
        term_a = 10 * torch.rand(self.batch_size) * sign
        term_b = 20 * torch.rand(self.batch_size) - 10
        term_y = term_a.unsqueeze(dim=1) * (
            torch.linspace(-16, 16, steps=32) ** 2
        ) + term_b.unsqueeze(dim=1)

        return {"input": term_y, "target": is_positive}


class ExampleDataModule(LightningDataModule):
    """
    An example data module for Pytorch Lightning.
    Args:
        batch_size: batch size
    """

    def __init__(self, batch_size: int = 64):
        super().__init__()
        self.batch_size = batch_size

    def prepare_data(self) -> None:
        pass

    def setup(self, stage: str):
        if stage == "fit":
            self.train_dataset = ExampleDataset(batch_size=self.batch_size)

        if stage == "val" or stage == "fit":
            self.val_dataset = ExampleDataset(batch_size=self.batch_size)

        if stage == "test":
            self.test_dataset = ExampleDataset(batch_size=self.batch_size)

        # log some "hyperparameters" that are for some reason
        # only accessible in the data module. An example could be
        # to store a hash of the loaded files for reproducibility.
        my_file_hash = "1234567890"
        self.trainer.logger.log_hyperparams(
            {"hyper_parameters": {"data_file_hash": my_file_hash}}
        )

    def train_dataloader(self):
        # batch_size = None signals to the DataLoader that
        # the Dataset is already batched
        return DataLoader(self.train_dataset, batch_size=None)

    def val_dataloader(self):
        return DataLoader(self.val_dataset, batch_size=None)

    def test_dataloader(self):
        return DataLoader(self.test_dataset, batch_size=None)


class ExampleNetwork(torch.nn.Module):
    """
    Example Network returning logits
    """

    def __init__(
        self,
        input_size=32,
        hidden_size=32,
        output_size=1,
    ):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(input_size, hidden_size),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_size, hidden_size),
            torch.nn.ReLU(),
            torch.nn.Linear(hidden_size, output_size),
        )

    def forward(self, batch):
        logits = self.net(batch["input"]).squeeze()
        return {"logits": logits, "target": batch["target"]}


class ExampleLightningModel(pl.LightningModule):
    def __init__(self):
        super().__init__()
        self.model = ExampleNetwork()
        self.criterion = torch.nn.BCEWithLogitsLoss()

    def forward(self, batch):
        return self.model(batch)

    def step(self, batch, subset):
        outputs = self(batch)
        loss = self.criterion(outputs["logits"], batch["target"].float())
        self.log(f"{subset}_loss", loss)
        return loss

    def training_step(self, batch, batch_idx):
        return self.step(batch, "train")

    def validation_step(self, batch, batch_idx):
        return self.step(batch, "val")

    def test_step(self, batch, batch_idx):
        return self.step(batch, "test")


def cli_main(args: ArgsType = None):
    LightningCLI(
        model_class=ExampleLightningModel,
        datamodule_class=ExampleDataModule,
        save_config_callback=LaminSaveConfigCallback,
        args=args,
    )


if __name__ == "__main__":
    cli_main()
