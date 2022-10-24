import sklearn.datasets
from lndb_setup import init, settings
from lndb_setup._settings_store import InstanceSettingsStore

import lamindb as ln
from lamindb.schema import core


def test_dynamic_settings():
    init(storage="test-instance-1", dbconfig="sqlite", schema="bionty")

    settings_store = InstanceSettingsStore(
        storage_root=str(settings.instance.storage_root),
        storage_region=settings.instance.storage_region,
        schema_modules=settings.instance.schema_modules,
        dbconfig=settings.instance._dbconfig,
    )

    pipeline = ln.db.add(core.pipeline(v="1", name="test-pipeline"))
    pipeline_run = ln.schema.core.pipeline_run(
        pipeline_id=pipeline.id, pipeline_v=pipeline.v, name="test-run"
    )

    ingest = ln.db.Ingest(dsource=pipeline_run)
    df = sklearn.datasets.load_iris(as_frame=True).frame
    ingest.add(df, name="test-dobject-1")
    ingest.commit()

    select_dobject_result = ln.db.select(core.dobject, name="test-dobject-1").all()
    assert len(select_dobject_result) == 1

    init(storage="test-instance-2", dbconfig="sqlite")

    select_dobject_result = ln.db.select(core.dobject, name="test-dobject-1").all()
    assert len(select_dobject_result) == 0

    select_dobject_result = ln.db.select(
        core.dobject, _settings_store=settings_store, name="test-dobject-1"
    ).all()
    assert len(select_dobject_result) == 1

    db_metadata = ln.schema._core.get_db_metadata_as_dict(settings_store)
    db_path = "/home/runner/work/lamindb/lamindb/test-instance-1/test-instance-1"
    assert db_metadata["key"] == db_path
