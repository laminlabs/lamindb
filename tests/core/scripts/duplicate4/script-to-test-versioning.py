import lamindb as ln
import pytest
from lamindb.core.exceptions import MissingContextUID

with pytest.raises(MissingContextUID) as error:
    ln.track()

assert (
    "there are already multiple transforms whose keys end with filename"
    in error.exconly()
)

ln.track("PC0eB2QPm0jW0000")
