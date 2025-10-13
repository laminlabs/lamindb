import lamindb as ln
import lightning as pl
import pytest
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


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


def test_callback_basic():
    artifact_key = "test/model.ckpt"

    callback = ln.integrations.lightning.Callback(
        path="test_checkpoint.ckpt", key=artifact_key
    )

    model = SimpleModel()
    train_data = DataLoader(
        TensorDataset(torch.randn(100, 10), torch.randn(100, 1)), batch_size=10
    )

    trainer = pl.Trainer(
        max_epochs=2, callbacks=[callback], enable_checkpointing=False, logger=False
    )
    trainer.fit(model, train_data)

    artifacts = ln.Artifact.filter(key=artifact_key).all()
    assert len(artifacts) == 2

    for af in artifacts:
        assert af.kind == "model"
        af.delete(permanent=True)


def test_callback_with_features():
    ln.Feature(name="train_loss", dtype="float").save()
    ln.Feature(name="custom_param", dtype="str").save()

    artifact_key = "test/model_features.ckpt"

    callback = ln.integrations.lightning.Callback(
        path="test_checkpoint_features.ckpt",
        key=artifact_key,
        features={"train_loss": None, "custom_param": "test_value"},
    )

    model = SimpleModel()
    train_data = DataLoader(
        TensorDataset(torch.randn(100, 10), torch.randn(100, 1)), batch_size=10
    )

    trainer = pl.Trainer(
        max_epochs=2, callbacks=[callback], enable_checkpointing=False, logger=False
    )
    trainer.fit(model, train_data)

    artifacts = ln.Artifact.filter(key=artifact_key).all()
    assert len(artifacts) == 2

    for af in artifacts:
        values = af.features.get_values()
        assert "train_loss" in values
        assert values["custom_param"] == "test_value"
        af.delete(permanent=True)


def test_callback_missing_features():
    artifact_key = "test/model_missing.ckpt"

    callback = ln.integrations.lightning.Callback(
        path="test_checkpoint_missing.ckpt",
        key=artifact_key,
        features={"nonexistent_feature": None},
    )

    model = SimpleModel()
    train_data = DataLoader(
        TensorDataset(torch.randn(100, 10), torch.randn(100, 1)), batch_size=10
    )

    trainer = pl.Trainer(
        max_epochs=1, callbacks=[callback], enable_checkpointing=False, logger=False
    )

    with pytest.raises(ValueError, match="Feature nonexistent_feature missing"):
        trainer.fit(model, train_data)
