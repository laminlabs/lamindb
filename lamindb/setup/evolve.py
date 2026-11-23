"""Make the setup backend available under the public dotted import path."""

import sys

from lamindb_setup import evolve

sys.modules[__name__] = evolve
