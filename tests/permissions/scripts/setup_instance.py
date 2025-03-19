import lamindb_setup as ln_setup
from hubmodule._setup import _install_db_module
from laminci.db import setup_local_test_postgres

pgurl = setup_local_test_postgres()

ln_setup.init(
    storage="./default_storage_permissions",
    name="lamindb-test-permissions",
    db=pgurl,
)

_install_db_module(pgurl, public_instance=False)
# can't add this app in the init because don't want t trigger the initial migration
# that conflicts with _install_db_module
ln_setup.settings.instance._schema_str = "hubmodule"
ln_setup.settings.instance._persist()
