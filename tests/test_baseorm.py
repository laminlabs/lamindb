import pytest

import lamindb as ln


def test_init_with_args():
    with pytest.raises(ValueError):
        ln.Tag("an arg")


def test_validate_required_fields():
    # tag has a required name
    with pytest.raises(TypeError):
        ln.Tag()
    # project has a required name
    with pytest.raises(TypeError):
        ln.Project(external_id="test")
