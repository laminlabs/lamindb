import subprocess
from pathlib import Path
from typing import Optional

from lamin_utils import logger
from lamindb_setup.core.hashing import hash_code

from ._settings import settings


def clone_git_repo(git_url: str) -> None:
    if Path(git_url).exists():
        logger.warning("git repo already exists")
        return None
    if not git_url.endswith(".git"):
        git_url += ".git"
    logger.important(f"cloning {git_url}")
    result = subprocess.run(
        f"git clone --depth 10 {git_url}",
        shell=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.decode())


def dir_from_repo_url(repo_url: Optional[str]) -> Optional[str]:
    if repo_url is not None:
        repo_dir = repo_url.split("/")[-1].replace(".git", "")
    return repo_dir


def get_git_commit_hash(
    blob_hash: str, cd_repo: Optional[str]
) -> subprocess.CompletedProcess:
    return subprocess.run(
        f"git log --find-object={blob_hash} --pretty=format:%H",
        shell=True,
        capture_output=True,
        cwd=cd_repo,
    )


def get_filepath_within_git_repo(
    commit_hash: str, blob_hash: str, cd_repo: Optional[str]
) -> str:
    result = subprocess.run(
        f"git ls-tree -r {commit_hash} | grep -E {blob_hash}",
        shell=True,
        capture_output=True,
        cwd=cd_repo,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"git ls-tree -r {commit_hash} | grep -E {blob_hash}\n"
            + result.stderr.decode()
        )
    if len(result.stdout.decode()) == 0:
        raise RuntimeError("Could not find filepath within git repo.")
    filepath = result.stdout.decode().split()[-1]
    return filepath


def get_transform_reference_from_git_repo(path: Path):
    blob_hash = hash_code(path).hexdigest()
    cd_repo = None
    result = get_git_commit_hash(blob_hash, cd_repo=None)
    commit_hash = result.stdout.decode()
    if commit_hash == "" or result.returncode == 1:
        cd_repo = dir_from_repo_url(settings.sync_git_repo)
        clone_git_repo(settings.sync_git_repo)
        result = get_git_commit_hash(blob_hash, cd_repo=cd_repo)
        commit_hash = result.stdout.decode()
        if commit_hash == "" or result.returncode == 1:
            raise RuntimeError(
                f"Did not find file in git repo\n{result.stderr.decode()}"
            )
    gitpath = get_filepath_within_git_repo(commit_hash, blob_hash, cd_repo)
    reference = f"{settings.sync_git_repo}/blob/{commit_hash}/{gitpath}"
    reference_type = "url"
    return reference, reference_type
