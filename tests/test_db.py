from pathlib import Path

import lamindb as lndb
from lamindb import setup


def test_create_to_load():
    storage = Path.home() / "mydata"
    setup.setup_from_cli(
        storage=storage,
        user="raspbear@gmx.de",
        secret="sd9lsxH6r5gjwKRSftTc9LNnIP5ua6Hhw8TVxl6l",
    )

    lndb.admin.db.insert.dobject("test_file", ".csv", interface_id="83jf")
    for entity in lndb.track.schema.entities:
        lndb.do.load(entity)

    (storage / "mydata.lndb").unlink()


if __name__ == "__main__":
    test_create_to_load()
