"""Universal IDs.

Base generators:

.. autosummary::
   :toctree: .

   base26
   base62
   base64

`uid` generators:

.. autosummary::
   :toctree: .

   base62_8
   base62_12
   base62_16
   base62_20


Collision probabilities
=======================

8 base62 characters (`62**8=2e+14`):

======= ===========
n       p_collision
======= ===========
100k    2e-05
1M      2e-03
======= ===========

12 base62 characters (`62**12=3e+21`):

======= ===========
n       p_collision
======= ===========
100M    2e-06
1B      2e-04
======= ===========

16 base62 characters (`62**16=5e+28`):

======= ===========
n       p_collision
======= ===========
1e12    7e-05
1e13    7e-03
======= ===========

20 base62 characters (`62**20=7e+35`) roughly matches UUID (`2**122=5e+36`):

======= ===========
n       p_collision
======= ===========
1e16    7e-05
1e17    7e-03
======= ===========

See `source <https://lamin.ai/laminlabs/lamindata/transform/t2xCdMB9v5wL>`__.

"""

import secrets
import string


def base64(n_char: int) -> str:
    """Random Base64 string."""
    alphabet = string.digits + string.ascii_letters.swapcase() + "_" + "-"
    uid = "".join(secrets.choice(alphabet) for i in range(n_char))
    return uid


def base62(n_char: int) -> str:
    """Random Base62 string."""
    alphabet = string.digits + string.ascii_letters.swapcase()
    uid = "".join(secrets.choice(alphabet) for i in range(n_char))
    return uid


def base26(n_char: int):
    """ASCII lowercase."""
    alphabet = string.ascii_lowercase
    uid = "".join(secrets.choice(alphabet) for i in range(n_char))
    return uid


def base62_4() -> str:
    return base62(4)


def base62_8() -> str:
    """Random Base62 string of length 8."""
    return base62(8)


def base62_12() -> str:
    """Random Base62 string of length 12."""
    return base62(12)


def base62_16() -> str:
    """Random Base62 string of length 16."""
    return base62(16)


def base62_20() -> str:
    """Random Base62 string of length 20."""
    return base62(20)


def base62_24() -> str:
    """Random Base62 string of length 24."""
    return base62(24)
