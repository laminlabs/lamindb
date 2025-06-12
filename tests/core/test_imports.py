import subprocess
import sys
import tempfile

import lamindb_setup.core.django as django_lamin
import pytest

TOPLEVEL_MODEL_CLASSES = [
    "Artifact",
    "Branch",
    "Collection",
    "Feature",
    "Person",
    "Project",
    "Record",
    "Reference",
    "Run",
    "Schema",
    "Sheet",
    "Space",
    "Storage",
    "Transform",
    "ULabel",
    "User",
    "FeatureSet",
    "Param",
]


@pytest.fixture(scope="module")
def connected_instance():
    if not django_lamin.IS_SETUP:
        pytest.xfail("Instance is not connected.")


@pytest.mark.usefixtures("connected_instance")
def test_ensure_exported_models():
    # ensure that we test all exported Django models in this test module
    import lamindb as ln
    from django.db.models import Model

    all_exported_django_models = {
        cls.__name__ for key, cls in vars(ln).items() if issubclass(cls, Model)
    }
    assert set(TOPLEVEL_MODEL_CLASSES) == all_exported_django_models


PRE_CONNECT_MODEL_TEMPLATE = """\
from lamindb import {model_class}  # import before connect
from lamindb import connect
connect("{instance_slug}")
{model_class}.df()
"""


@pytest.fixture(scope="function")
def disable_auto_connect():
    """Disable auto-connect for the duration of the test."""
    from lamindb_setup import settings

    original_auto_connect = settings.auto_connect
    settings.auto_connect = False
    try:
        yield
    finally:
        settings.auto_connect = original_auto_connect


@pytest.mark.parametrize("clsname", TOPLEVEL_MODEL_CLASSES)
@pytest.mark.usefixtures("disable_auto_connect")
def test_ensure_preconnect_model_import(clsname):
    from lamindb import settings

    with tempfile.NamedTemporaryFile(suffix=".py") as f:
        f.write(
            PRE_CONNECT_MODEL_TEMPLATE.format(
                model_class=clsname,
                instance_slug=settings.instance.slug,
            )
        )

        proc = subprocess.run([sys.executable, f.name])
    assert proc.returncode == 0
