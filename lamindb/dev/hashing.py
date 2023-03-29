import base64
import hashlib
from pathlib import Path
from typing import Set


def to_b64_str(bstr: bytes):
    b64 = base64.urlsafe_b64encode(bstr).decode().strip("=")
    return b64


# a lot to read about this: lamin-notes/2022/hashing
def hash_set(s: Set[str]) -> str:
    bstr = ":".join(sorted(s)).encode("utf-8")
    # as we're truncating at 20 b64, we choose md5 over sha512
    return to_b64_str(hashlib.md5(bstr).digest())[:20]


def hash_file(path: Path) -> str:
    # based on https://stackoverflow.com/questions/3431825/generating-an-md5-hash-of-a-file  # noqa
    hash_md5 = hashlib.md5()
    with open(path, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash_md5.update(chunk)
    return to_b64_str(hash_md5.digest())
