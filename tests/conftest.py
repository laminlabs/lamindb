import shutil
from pathlib import Path

import pytest


@pytest.fixture(scope="function")
def clean_soma_files(request):
    path = request.param if hasattr(request, "param") else "small_dataset.tiledbsoma"
    if Path(path).exists():
        shutil.rmtree(path)

    yield path

    if Path(path).exists():
        shutil.rmtree(path)
