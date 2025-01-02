"""IDs.

Base generators:

.. autosummary::
   :toctree: .

   base26
   base62
   base64

8 base62 characters:

======= ===========
n       p_collision
======= ===========
100k    2e-05
1M      2e-03
======= ===========

12 base62 characters:

======= ===========
n       p_collision
======= ===========
100M    2e-06
1B      2e-04
======= ===========

20 base62 characters (62**20=7e+35) roughly matches UUID (2*122=5e+36):

======= ===========
n       p_collision
======= ===========
3e15    1e-6
======= ===========

"""

import secrets
import string


def base64(n_char: int) -> str:
    """Random Base64 string."""
    alphabet = string.digits + string.ascii_letters.swapcase() + "_" + "-"
    id = "".join(secrets.choice(alphabet) for i in range(n_char))
    return id


def base62(n_char: int) -> str:
    """Random Base62 string."""
    alphabet = string.digits + string.ascii_letters.swapcase()
    id = "".join(secrets.choice(alphabet) for i in range(n_char))
    return id


# the following cannot be serialized by Django
# class Base62:
#     def __init__(self, n_char: int):
#         self.n_char = n_char

#     def __call__(self):
#         return base62(self.n_char)


def base26(n_char: int):
    """ASCII lowercase."""
    alphabet = string.ascii_lowercase
    id = "".join(secrets.choice(alphabet) for i in range(n_char))
    return id


def base62_4() -> str:
    return base62(4)


def base62_8() -> str:
    return base62(8)


def base62_12() -> str:
    return base62(12)


def base62_14() -> str:
    return base62(14)


def base62_16() -> str:
    return base62(16)


def base62_18() -> str:
    return base62(18)


def base62_20() -> str:
    return base62(20)


def base62_24() -> str:
    return base62(24)
