from pathlib import Path

import lamindb as lndb
from lamindb import setup


def test_create_to_load():
    setup.setup(
        storage=Path.home() / "data",
        user="falexwolf",
    )

    lndb.admin.db.insert.file("test_file.csv", interface_id="83jf")
    for entity in lndb.track.schema.entities:
        lndb.do.load(entity)


if __name__ == "__main__":
    test_create_to_load()
