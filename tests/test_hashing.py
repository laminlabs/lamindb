import base64
from pathlib import Path

from lamindb.dev.hashing import hash_file, to_b64_str


def test_compute_hash():
    dataset = [
        ("file_1.txt", "a", "DMF1ucDxtqgxw5niaXcmYQ"),
    ]
    for path, content, hash in dataset:
        filepath = Path(path)
        with open(path, "w") as file:
            file.write(content)
        assert hash == hash_file(filepath)
        filepath.unlink()


def test_base64():
    # the following can be commented out over time
    mytest = "test".encode()
    b64_str = to_b64_str(mytest)
    b64_str_padded = f"{b64_str}=="
    assert base64.urlsafe_b64decode(b64_str_padded.encode()).hex() == mytest.hex()
