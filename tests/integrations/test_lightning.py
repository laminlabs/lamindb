from pathlib import Path
from typing import Callable, Generator

import lamindb as ln
import lightning as pl
import pytest
import torch
from lamindb.integrations import lightning as ll
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


@pytest.fixture(autouse=True)
def cleanup_checkpoints() -> Generator[None, None, None]:
    yield
    for path in Path().glob("test_checkpoint*.ckpt"):
        path.unlink(missing_ok=True)


@pytest.fixture
def simple_model() -> pl.LightningModule:
    class SimpleModel(pl.LightningModule):
        def __init__(self):
            super().__init__()
            self.layer = nn.Linear(10, 1)

        def forward(self, x):
            return self.layer(x)

        def training_step(self, batch, batch_idx):
            x, y = batch
            loss = nn.functional.mse_loss(self(x), y)
            self.log("train_loss", loss)
            return loss

        def configure_optimizers(self):
            return torch.optim.Adam(self.parameters())

    model = SimpleModel()
    return model


@pytest.fixture
def torch_train_data_dataloader() -> DataLoader:
    train_data = DataLoader(
        TensorDataset(torch.randn(100, 10), torch.randn(100, 1)), batch_size=10
    )
    return train_data


def test_callback_basic(
    cleanup_checkpoints: Callable,
    torch_train_data_dataloader: DataLoader,
    simple_model: pl.LightningModule,
):
    """Callback should create artifacts for each training epoch."""
    artifact_key = "test/model.ckpt"

    callback = ll.Callback(path="test_checkpoint.ckpt", key=artifact_key)

    trainer = pl.Trainer(
        max_epochs=2, callbacks=[callback], enable_checkpointing=False, logger=False
    )
    trainer.fit(simple_model, torch_train_data_dataloader)

    artifacts = ln.Artifact.filter(key=artifact_key)
    assert len(artifacts) == 2

    for af in artifacts:
        assert af.kind == "model"
        af.delete(permanent=True)


def test_callback_with_features(
    cleanup_checkpoints: Callable,
    torch_train_data_dataloader: DataLoader,
    simple_model: pl.LightningModule,
):
    """Callback should annotate artifacts with feature values."""
    train_loss = ln.Feature(name="train_loss", dtype="float").save()
    custom_param = ln.Feature(name="custom_param", dtype="str").save()

    artifact_key = "test/model_features.ckpt"

    callback = ll.Callback(
        path="test_checkpoint_features.ckpt",
        key=artifact_key,
        features={"train_loss": None, "custom_param": "test_value"},
    )

    trainer = pl.Trainer(
        max_epochs=2, callbacks=[callback], enable_checkpointing=False, logger=False
    )
    trainer.fit(simple_model, torch_train_data_dataloader)

    artifacts = ln.Artifact.filter(key=artifact_key)
    assert len(artifacts) == 2

    for af in artifacts:
        values = af.features.get_values()
        assert "train_loss" in values
        assert values["custom_param"] == "test_value"
        af.delete(permanent=True)

    train_loss.delete(permanent=True)
    custom_param.delete(permanent=True)


def test_callback_missing_features(
    cleanup_checkpoints: Callable,
    torch_train_data_dataloader: DataLoader,
    simple_model: pl.LightningModule,
):
    """Callback should raise an error when specified features do not exist."""
    artifact_key = "test/model_missing.ckpt"

    callback = ll.Callback(
        path="test_checkpoint_missing.ckpt",
        key=artifact_key,
        features={"nonexistent_feature": None},
    )

    trainer = pl.Trainer(
        max_epochs=1, callbacks=[callback], enable_checkpointing=False, logger=False
    )

    with pytest.raises(ValueError) as e:
        trainer.fit(simple_model, torch_train_data_dataloader)
    assert "Feature nonexistent_feature missing" in str(e.value)
