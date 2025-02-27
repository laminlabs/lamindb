import lamindb_setup as ln_setup
from hubmodule._setup import _install_db_module
from laminci.db import setup_local_test_postgres

pgurl = setup_local_test_postgres()

ln_setup.init(
    storage="./default_storage_permissions",
    name="lamindb-test-permissions",
    db=pgurl,
)

_install_db_module(pgurl)

ln_setup.settings.instance._schema_str = "hubmodule"
ln_setup.settings.instance._persist()
