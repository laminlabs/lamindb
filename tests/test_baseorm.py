import pytest

import lamindb as ln


def test_init_with_args():
    with pytest.raises(ValueError):
        ln.Tag()


def test_validate_required_fields():
    # project has a required name
    with pytest.raises(TypeError):
        ln.Project(external_id="test")
