import subprocess
from pathlib import Path
from typing import Optional

from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.hashing import hash_code

from ._settings import sanitize_git_repo_url, settings


def get_git_repo_from_remote() -> Path:
    repo_url = settings.sync_git_repo
    repo_dir = setup_settings.storage.cache_dir / repo_url.split("/")[-1]
    if repo_dir.exists():
        logger.warning(f"git repo {repo_dir} already exists locally")
        return None
    logger.important(f"cloning {repo_url} into {repo_dir}")
    result = subprocess.run(
        f"git clone --depth 10 {repo_url}.git",
        shell=True,
        capture_output=True,
        cwd=setup_settings.storage.cache_dir,
    )
    if result.returncode != 0 or not repo_dir.exists():
        raise RuntimeError(result.stderr.decode())
    return repo_dir


def check_remote_git_url_matches_setting():
    result = subprocess.run(
        "git config --get remote.origin.url",
        shell=True,
        capture_output=True,
    )
    remote_url = sanitize_git_repo_url(result.stdout.decode().strip())
    assert remote_url == settings.sync_git_repo


def get_git_commit_hash(
    blob_hash: str, repo_dir: Optional[Path] = None
) -> Optional[str]:
    command = f"git log --find-object={blob_hash} --pretty=format:%H"
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        cwd=repo_dir,
    )
    commit_hash = result.stdout.decode()
    if commit_hash == "" or result.returncode == 1:
        return None
    else:
        assert len(commit_hash) == 40
        return commit_hash


def get_filepath_within_git_repo(
    commit_hash: str, blob_hash: str, repo_dir: Optional[Path]
) -> str:
    # repo_dir might not point to the root of the
    # the git repository because git log --find-object works
    # from anywhere in the repo, hence, let's get the root
    repo_root = (
        subprocess.run(
            "git rev-parse --show-toplevel",
            shell=True,
            capture_output=True,
            cwd=repo_dir,
        )
        .stdout.decode()
        .strip()
    )
    command = f"git ls-tree -r {commit_hash} | grep -E {blob_hash}"
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        cwd=repo_root,
    )
    if result.returncode != 0 and result.stderr.decode() != "":
        raise RuntimeError(f"{command}\n{result.stderr.decode()}")
    if len(result.stdout.decode()) == 0:
        raise RuntimeError(
            f"Could not find path in git repo {settings.sync_git_repo} running:\n{command}"
            f"\nin local clone: {repo_root}"
        )
    filepath = result.stdout.decode().split()[-1]
    return filepath


def get_transform_reference_from_git_repo(path: Path) -> str:
    blob_hash = hash_code(path).hexdigest()
    commit_hash = get_git_commit_hash(blob_hash)
    if commit_hash is not None:
        logger.warning("found script in local repository")
        check_remote_git_url_matches_setting()
        repo_dir = None
    else:
        repo_dir = get_git_repo_from_remote()
        commit_hash = get_git_commit_hash(blob_hash, repo_dir=repo_dir)
        if commit_hash is None:
            raise RuntimeError(
                f"Did not find file hash for {path} in git repo {repo_dir}"
            )
    gitpath = get_filepath_within_git_repo(commit_hash, blob_hash, repo_dir)
    reference = f"{settings.sync_git_repo}/blob/{commit_hash}/{gitpath}"
    return reference
