from pathlib import Path

import lamindb as lndb
from lamindb._setup import _setup


def test_create_to_load():
    storage = Path.home() / "mydata"
    _setup.log_in_user(
        email="raspbear@gmx.de",
        secret="O6A5bqHzhERTMpGHXxEdbsTrEIt4a3Dy3AWdrWQR",
    )
    _setup.setup_instance(storage=storage)
    lndb.admin.db.insert.dobject(
        name="test_file",
        file_suffix=".csv",
        interface_id="83jf",
        interface_v="1",
        interface_name="test",
        interface_type="other",
    )
    for entity in lndb.track.schema.entities:
        lndb.do.load(entity)

    (storage / "mydata.lndb").unlink()


if __name__ == "__main__":
    test_create_to_load()
