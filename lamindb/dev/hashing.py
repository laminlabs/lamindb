"""Hashing.

.. autosummary::
   :toctree: .

   hash_set
   hash_file

"""

import base64
import hashlib
from typing import Set, Tuple


def to_b64_str(bstr: bytes):
    b64 = base64.urlsafe_b64encode(bstr).decode().strip("=")
    return b64


def b16_to_b64(s: str):
    return to_b64_str(base64.b16decode(s.strip('"'), casefold=True))


# a lot to read about this: lamin-notes/2022/hashing
def hash_set(s: Set[str]) -> str:
    bstr = ":".join(sorted(s)).encode("utf-8")
    # as we're truncating at 20 b64, we choose md5 over sha512
    return to_b64_str(hashlib.md5(bstr).digest())[:20]


def hash_file(file_path, chunk_size=50 * 1024 * 1024) -> Tuple[str, str]:
    chunks = []
    with open(file_path, "rb") as fp:
        # read first chunk
        chunks = [fp.read(chunk_size)]
        # try reading the 2nd chunk
        data = fp.read(chunk_size)
        if data:
            # go to end of file
            fp.seek(-chunk_size, 2)
            # read last chunk
            data = fp.read(chunk_size)
            chunks.append(data)
    if len(chunks) == 1:
        digest = hashlib.md5(chunks[0]).digest()
        hash_type = "md5"
    else:
        digests = b"".join(hashlib.sha1(chunk).digest() for chunk in chunks)
        digest = hashlib.sha1(digests).digest()
        hash_type = "sha1-fl"  # sha1 first last chunk
    return to_b64_str(digest)[:22], hash_type
