import base64
from pathlib import Path

from lamindb._record import hash_file, to_b64_str
from lamindb.dev._core import get_name_suffix_from_filepath


def test_get_name_suffix_from_filepath():
    # based on https://stackoverflow.com/questions/31890341/clean-way-to-get-the-true-stem-of-a-path-object  # noqa
    dataset = [
        ("a", "a", ""),
        ("a.txt", "a", ".txt"),
        ("archive.tar.gz", "archive", ".tar.gz"),
        ("directory/file", "file", ""),
        ("d.x.y.z/f.a.b.c", "f", ".a.b.c"),
        ("logs/date.log.txt", "date", ".txt"),
    ]
    for path, name, suffix in dataset:
        filepath = Path(path)
        assert name, suffix == get_name_suffix_from_filepath(filepath)


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
