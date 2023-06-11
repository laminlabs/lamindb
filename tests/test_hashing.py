import base64
from pathlib import Path

from lamindb.dev.hashing import hash_file, to_b64_str


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
    # the following can be commented out over time
    mytest = "test".encode()
    b64_str = to_b64_str(mytest)
    b64_str_padded = f"{b64_str}=="
    assert base64.urlsafe_b64decode(b64_str_padded.encode()).hex() == mytest.hex()
