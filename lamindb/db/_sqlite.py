import random
import string

from sqlalchemy import Column, MetaData, String, Table, create_engine


def uid(n_char=6):
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


class DB:
    """Interact with local database & storage.

    We work with base62 primary keys:

    ====== =========
    len_id n_entries
    ====== =========
    1      >6e+01
    2      >4e+03
    3      >2e+05
    4      >1e+07
    5      >9e+08
    6      >6e+10
    7      >4e+12
    8      >2e+14
    9      >1e+16
    12     >3e+21 (nbproject uid)
    20     >7e+35 (~UUID)
    ====== =========
    """

    def __init__(self, seed: int = 0):
        from lamindb._configuration import storage_root

        random.seed(seed)
        # fixed location of database
        self._database_file = storage_root + "lamin.db"

    def create(self, seed: int = 0):
        """Create database with initial schema."""
        # use the schema we just migrated to SQL and add a primary key
        metadata = MetaData()

        def _uid_file():
            return uid(n_char=20)

        Table(
            "file",
            metadata,
            Column("id", String, primary_key=True, default=_uid_file),
            Column("name", String),
            Column("source", String),  # can be anything, e.g. a notebook
        )

        Table(
            "notebook",
            metadata,
            Column("id", String, primary_key=True),  # this is an nbproject uid
            Column("source", String),  # can be anything, e.g. a URL
        )

        engine = create_engine(f"sqlite:///{self._database_file}", future=True)
        metadata.create_all(bind=engine)
