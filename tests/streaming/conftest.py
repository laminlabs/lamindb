import pytest


@pytest.fixture(scope="session", autouse=True)
def streaming_instance(tmp_path_factory):
    root = tmp_path_factory.mktemp("streaming-instance")
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("LAMIN_SETTINGS_DIR", str(root / "settings"))
        patch.setenv("LAMIN_CACHE_DIR", str(root / "cache"))
        patch.setenv("LAMIN_TESTING", "true")
        import lamindb as ln
        import lamindb_setup

        lamindb_setup._TESTING = True
        ln.setup.init(storage=root / "storage", name="streaming-prototype-tests")
        ln.settings.creation.artifact_silence_missing_run_warning = True
        yield
        ln.setup.settings.dev_dir = None
