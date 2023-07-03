import base64
from pathlib import Path

from lamindb.dev.hashing import b16_to_b64, hash_file, to_b64_str


def test_compute_hash():
    files = [
        # filepath, content, hash
        ("a", "DMF1ucDxtqgxw5niaXcmYQ", None),
        ("b", "kutf_uauL-w61xx3dTFXjw", None),
        ("abc", "kAFQmDzST7DWlj99KOF_cg", None),
        ("a", "DMF1ucDxtqgxw5niaXcmYQ", 1),
        ("b", "kutf_uauL-w61xx3dTFXjw", 1),
        # the last case here triggers multi-chunk compute with sha1
        ("abc", "od1ZbYc380Hkre4BIjariX", 1),
    ]
    for content, hash, chunk_size in files:
        filepath = Path("file_1.txt")
        filepath.write_text(content)
        assert hash == hash_file(filepath, chunk_size=chunk_size)
        filepath.unlink()


def test_base64():
    mytest = "test".encode()
    b64_str = to_b64_str(mytest)
    b64_str_padded = f"{b64_str}=="
    assert base64.urlsafe_b64decode(b64_str_padded.encode()).hex() == mytest.hex()


def test_b16_to_b64():
    assert b16_to_b64("9b89c8c1acf79dba5b5341d1fff9806f") == "m4nIwaz3nbpbU0HR__mAbw"
