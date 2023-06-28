import base64
from pathlib import Path

from lamindb.dev.hashing import b16_to_b64, hash_file, to_b64_str


def test_compute_hash():
    files = [
        # filepath, content, hash
        ("file_1.txt", "a", "DMF1ucDxtqgxw5niaXcmYQ"),
        ("file_1.txt", "b", "kutf_uauL-w61xx3dTFXjw"),
    ]
    for path, content, hash in files:
        filepath = Path(path)
        filepath.write_text(content)
        assert hash == hash_file(filepath)
        filepath.unlink()


def test_base64():
    mytest = "test".encode()
    b64_str = to_b64_str(mytest)
    b64_str_padded = f"{b64_str}=="
    assert base64.urlsafe_b64decode(b64_str_padded.encode()).hex() == mytest.hex()


def test_b16_to_b64():
    assert b16_to_b64("9b89c8c1acf79dba5b5341d1fff9806f") == "m4nIwaz3nbpbU0HR__mAbw"
