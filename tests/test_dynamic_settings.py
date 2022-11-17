import sklearn.datasets
from lndb_setup import init, settings
from lndb_setup._settings_store import InstanceSettingsStore

import lamindb as ln


def test_dynamic_settings():
    init(storage="test-instance-1", dbconfig="sqlite", schema="bionty")

    settings_store = InstanceSettingsStore(
        storage_root=str(settings.instance.storage_root),
        storage_region=settings.instance.storage_region,
        schema_modules=settings.instance.schema_modules,
        dbconfig=settings.instance._dbconfig,
    )

    pipeline = ln.add(ln.schema.Pipeline(v="1", name="test-pipeline"))
    run = ln.schema.Run(pipeline_id=pipeline.id, pipeline_v=pipeline.v, name="test-run")

    df = sklearn.datasets.load_iris(as_frame=True).frame
    dobject = ln.record(df, name="test-dobject-1", run=run)
    ln.add(dobject)

    select_dobject_result = ln.select(ln.schema.DObject, name="test-dobject-1").all()
    assert len(select_dobject_result) == 1

    init(storage="test-instance-2", dbconfig="sqlite")

    select_dobject_result = ln.select(ln.schema.DObject, name="test-dobject-1").all()
    assert len(select_dobject_result) == 0

    select_dobject_result = ln.select(
        ln.schema.DObject, _settings_store=settings_store, name="test-dobject-1"
    ).all()
    assert len(select_dobject_result) == 1

    db_metadata = ln.schema._core.get_db_metadata_as_dict(settings_store)
    db_substring = "lamindb/test-instance-1/test-instance-1.lndb"
    assert db_substring in str(db_metadata["key"])
