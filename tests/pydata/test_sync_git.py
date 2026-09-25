import subprocess

import pytest
from lamindb.core._sync_git import get_filepath_within_git_repo


def _git(*args, cwd):
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.mark.parametrize(
    "relpath", ["my analysis/run qc.py", "analyse/qualité.py", "scripts/run.py"]
)
def test_get_filepath_within_git_repo_returns_the_full_path(tmp_path, relpath):
    _git("init", "-q", cwd=tmp_path)
    script = tmp_path / relpath
    script.parent.mkdir(parents=True)
    script.write_text("print('hello')\n")
    _git("add", ".", cwd=tmp_path)
    _git(
        "-c",
        "user.name=test",
        "-c",
        "user.email=test@example.com",
        "commit",
        "-q",
        "-m",
        "add script",
        cwd=tmp_path,
    )
    commit_hash = _git("rev-parse", "HEAD", cwd=tmp_path)
    blob_hash = _git("hash-object", str(script), cwd=tmp_path)

    filepath = get_filepath_within_git_repo(commit_hash, blob_hash, tmp_path)

    assert filepath == relpath
