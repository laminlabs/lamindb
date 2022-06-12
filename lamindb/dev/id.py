import random
import string


def id(n_char: int = 6) -> str:
    base62 = string.digits + string.ascii_letters.swapcase()
    id = "".join(random.choice(base62) for i in range(n_char))
    return id


def id_file() -> str:
    return id(n_char=20)


def id_user() -> str:
    return id(n_char=3)


def id_track() -> str:
    return id(n_char=24)
