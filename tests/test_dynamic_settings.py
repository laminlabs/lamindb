from lndb_setup import init, settings
from lndb_setup._settings_store import InstanceSettingsStore

import lamindb as ln
from lamindb.schema import core


def test_dynamic_settings():
    settings_store = InstanceSettingsStore(
        storage_root=str(settings.instance.storage_root),
        storage_region=settings.instance.storage_region,
        schema_modules=settings.instance.schema_modules,
        dbconfig=settings.instance._dbconfig,
    )
    init(storage="another-instance", dbconfig="sqlite", schema="bionty")

    select_dobject_result = ln.db.select(
        core.dobject, _settings_store=settings_store, name="iris"
    ).all()
    assert len(select_dobject_result) == 1

    db_metadata = ln.schema._core.get_db_metadata_as_dict(settings_store)
    assert db_metadata["key"] == "mydata-test-db"
