import shutil
from pathlib import Path
from typing import Generator

import lamindb as ln
import lightning as pl
import pytest
import torch
from lamindb.integrations import lightning as ll
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


@pytest.fixture(autouse=True)
def cleanup_checkpoints() -> Generator[None, None, None]:
    """Clean up checkpoint files and directories after each test."""
    yield
    checkpoints_dir = Path("checkpoints")
    if checkpoints_dir.exists():
        shutil.rmtree(checkpoints_dir)


@pytest.fixture(autouse=True, scope="session")
def cleanup_test_dir() -> Generator[None, None, None]:
    """Clean up test directory after all tests."""
    yield
    for dirname in ("lightning_checkpoints", "test_lightning"):
        dirpath = Path(dirname)
        if dirpath.exists():
            shutil.rmtree(dirpath)


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

    return SimpleModel()


@pytest.fixture
def dataloader() -> DataLoader:
    return DataLoader(
        TensorDataset(torch.randn(100, 10), torch.randn(100, 1)), batch_size=10
    )


@pytest.fixture
def dirpath(request: pytest.FixtureRequest) -> Generator[str, None, None]:
    prefix = f"lightning_checkpoints/{request.node.name}/"
    resolved = str(Path(prefix).resolve()) + "/"

    yield prefix

    for af in ln.Artifact.filter(key__startswith=resolved):
        af.delete(permanent=True, storage=True)
    dirpath_path = Path(prefix)
    if dirpath_path.exists():
        shutil.rmtree(dirpath_path)


@pytest.fixture(scope="session")
def lightning_features() -> Generator[None, None, None]:
    """Create lightning features."""
    ll.save_lightning_features()

    yield

    if lightning_type := ln.Feature.filter(name="lamindb.lightning").one_or_none():
        for feat in ln.Feature.filter(type=lightning_type):
            for af in ln.Artifact.filter(schemas__features=feat):
                af.delete(permanent=True, storage=True)
            feat.delete(permanent=True)


def test_checkpoint_basic(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
):
    """Checkpoint should create artifacts with semantic paths."""
    callback = ll.Checkpoint(dirpath=dirpath, monitor="train_loss")
    trainer = pl.Trainer(
        max_epochs=2,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    assert len(artifacts) >= 1
    for af in artifacts:
        assert af.kind == "model"
        assert af.key.startswith(resolved)


def test_checkpoint_with_features(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
):
    """Checkpoint should annotate artifacts with feature values."""
    ln.Feature(name="train_loss", dtype=float).save()
    ln.Feature(name="custom_param", dtype=str).save()

    ln.track()

    callback = ll.Checkpoint(
        dirpath=dirpath,
        features={
            "artifact": {"train_loss": None},
            "run": {"custom_param": "test_value"},
        },
        monitor="train_loss",
    )
    trainer = pl.Trainer(
        max_epochs=2,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    assert len(artifacts) >= 1
    for af in artifacts:
        values = af.features.get_values()
        assert "train_loss" in values

    assert ln.context.run.features.get_values()["custom_param"] == "test_value"

    ln.finish()


def test_checkpoint_missing_features(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
):
    """Checkpoint should raise an error when specified features do not exist."""
    callback = ll.Checkpoint(
        dirpath=dirpath,
        features={"artifact": {"nonexistent_feature": None}},
        monitor="train_loss",
    )
    trainer = pl.Trainer(
        max_epochs=1,
        callbacks=[callback],
        logger=False,
    )

    with pytest.raises(ValueError, match="Feature nonexistent_feature missing"):
        trainer.fit(simple_model, dataloader)


def test_checkpoint_auto_features(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
    lightning_features: None,
):
    """Checkpoint should auto-track lightning features if they exist."""
    callback = ll.Checkpoint(
        dirpath=dirpath,
        monitor="train_loss",
        save_top_k=2,
    )
    trainer = pl.Trainer(
        max_epochs=3,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    assert len(artifacts) >= 1

    for af in artifacts:
        values = af.features.get_values()
        assert "is_best_model" in values
        assert "score" in values
        assert "model_rank" in values


def test_checkpoint_best_model_tracking(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
    lightning_features: None,
):
    """Only one checkpoint should be marked as best model."""
    callback = ll.Checkpoint(
        dirpath=dirpath,
        monitor="train_loss",
        save_top_k=3,
        mode="min",
    )
    trainer = pl.Trainer(
        max_epochs=3,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    best_count = sum(
        1 for af in artifacts if af.features.get_values().get("is_best_model") is True
    )
    assert best_count == 1


def test_checkpoint_model_rank(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
    lightning_features: None,
):
    """Checkpoints should have correct model_rank (0 = best)."""
    callback = ll.Checkpoint(
        dirpath=dirpath,
        monitor="train_loss",
        save_top_k=3,
        mode="min",
    )
    trainer = pl.Trainer(
        max_epochs=3,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    ranks = [af.features.get_values().get("model_rank") for af in artifacts]
    assert 0 in ranks  # best model has rank 0


def test_checkpoint_semantic_paths(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
    lightning_features: None,
):
    """Checkpoints should have semantic keys derived from dirpath."""
    callback = ll.Checkpoint(
        dirpath=dirpath,
        monitor="train_loss",
        save_top_k=3,
    )
    trainer = pl.Trainer(
        max_epochs=3,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    resolved = callback.dirpath.rstrip("/") + "/"
    artifacts = ln.Artifact.filter(key__startswith=resolved)
    assert len(artifacts) >= 1

    for af in artifacts:
        assert af.key.startswith(resolved)
        values = af.features.get_values()
        assert "is_best_model" in values
        assert "score" in values


def test_callback_deprecated(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    tmp_path: Path,
):
    """Deprecated Callback should still work."""
    key = f"test/legacy/{tmp_path.name}/model.ckpt"
    path = tmp_path / "model.ckpt"

    with pytest.warns(DeprecationWarning, match="use ll.Checkpoint instead"):
        callback = ll.Callback(path=path, key=key)

    trainer = pl.Trainer(
        max_epochs=1,
        callbacks=[callback],
        logger=False,
    )
    trainer.fit(simple_model, dataloader)

    artifacts = ln.Artifact.filter(key=key)
    assert len(artifacts) >= 1
    assert artifacts[0].kind == "model"

    # cleanup
    for af in artifacts:
        af.delete(permanent=True, storage=True)


@pytest.mark.parametrize(
    "overwrite_versions,should_raise", [(False, True), (True, False)]
)
def test_checkpoint_overwrite(
    simple_model: pl.LightningModule,
    dataloader: DataLoader,
    dirpath: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    overwrite_versions: bool,
    should_raise: bool,
):
    """Checkpoint should raise when artifact exists unless overwrite_versions=True."""
    dummy = tmp_path / "dummy.ckpt"
    dummy.write_bytes(b"dummy")
    fixed_key = f"{dirpath.rstrip('/')}/fixed.ckpt"
    ln.Artifact(dummy, key=fixed_key).save()

    callback = ll.Checkpoint(dirpath=dirpath, overwrite_versions=overwrite_versions)
    monkeypatch.setattr(callback, "_get_artifact_key", lambda _: fixed_key)
    trainer = pl.Trainer(max_epochs=1, callbacks=[callback], logger=False)

    if should_raise:
        with pytest.raises(ValueError, match="already exists"):
            trainer.fit(simple_model, dataloader)

    else:
        trainer.fit(simple_model, dataloader)

    for af in ln.Artifact.filter(key=fixed_key):
        af.delete(permanent=True, storage=True)


def test_checkpoint_invalid_feature_keys(dirpath: str):
    """Checkpoint should raise on invalid feature keys."""
    with pytest.raises(ValueError, match="Invalid feature keys"):
        ll.Checkpoint(
            dirpath=dirpath,
            features={"invalid_key": {"foo": "bar"}},  # type: ignore
        )
