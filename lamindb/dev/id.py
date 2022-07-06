import random
import string


def id(n_char: int) -> str:
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def id_dobject() -> str:
    """IDs for dobject.

    Fix this table!

    ======= ===========
    n_users p_collision
    ======= ===========
    100k    2e-09
    1M      2e-07
    """
    return id(n_char=20)


def id_user() -> str:
    """User ID with 8 base62 characters.

    Consistent with 1M users producing 1k notebooks.
    Safe for 100k users producing 10k notebooks.

    Allows >2e14 users.

    This is one of 2 IDs that are centralized.

    Collision probability in decentralized system is:

    ======= ===========
    n_users p_collision
    ======= ===========
    100k    2e-05
    1M      2e-03
    """
    return id(n_char=8)


def id_track() -> str:
    return id(n_char=24)


def id_secret() -> str:
    return id(n_char=40)
