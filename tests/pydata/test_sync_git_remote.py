import subprocess

import lamindb as ln
import pytest
from lamindb.core._sync_git import check_local_git_repo


@pytest.fixture
def sync_git_repo():
    ln.settings.sync_git_repo = "https://github.com/laminlabs/lamindb"
    yield ln.settings.sync_git_repo
    ln.settings._sync_git_repo = None


@pytest.mark.parametrize(
    "remote_url,expected",
    [
        ("https://github.com/laminlabs/lamindb.git", True),
        ("git@github.com:laminlabs/lamindb.git", True),
        ("ssh://git@github.com/laminlabs/lamindb.git", True),
        ("git@github.com:laminlabs/lamin-cli.git", False),
        ("/srv/git/lamindb.git", False),
    ],
)
def test_check_local_git_repo_accepts_ssh_remotes(
    tmp_path, monkeypatch, sync_git_repo, remote_url, expected
):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(
        ["git", "remote", "add", "origin", remote_url], cwd=tmp_path, check=True
    )
    monkeypatch.chdir(tmp_path)

    assert check_local_git_repo() is expected
