import os


def test_nbconvert():
    exit_code = os.system(  # noqa: S605
        "jupyter nbconvert --to notebook --inplace --execute ./tests/core/notebooks/load_schema.ipynb"
    )
    assert exit_code == 0
