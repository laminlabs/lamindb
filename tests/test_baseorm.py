import pytest

import lamindb as ln


def test_init_with_args():
    with pytest.raises(ValueError):
        ln.Tag()
