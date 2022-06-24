import random
import string


def id(n_char: int) -> str:
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def id_file() -> str:
    return id(n_char=20)


def id_user() -> str:
    """User ID.

    Consistent with 1M users producing 1k notebooks.
    Safe for 100k users producing 10k notebooks.

    Allows >2e14 users.

    Collision probability in decentralized system is:

    ======= ===========
    n_users p_collision
    ======= ===========
    10k     2e-07
    100k    2e-05
    1M      2e-03
    """
    return id(n_char=8)


def id_track() -> str:
    return id(n_char=24)
