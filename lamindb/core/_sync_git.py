import subprocess
from pathlib import Path
from typing import Optional

from lamin_utils import logger
from lamindb_setup.core.hashing import hash_code

from ._settings import settings


def clone_git_repo(git_url: str) -> None:
    repo_dir = dir_from_repo_url(git_url)
    if Path(repo_dir).exists():
        logger.warning(f"git repo {repo_dir} already exists locally")
        return None
    if not git_url.endswith(".git"):
        git_url += ".git"
    logger.important(f"cloning {git_url}")
    result = subprocess.run(
        f"git clone --depth 10 {git_url}",
        shell=True,
        capture_output=True,
    )
    if result.returncode != 0 or not Path(repo_dir).exists():
        raise RuntimeError(result.stderr.decode())


def dir_from_repo_url(repo_url: Optional[str]) -> Optional[str]:
    if repo_url is not None:
        return repo_url.split("/")[-1].replace(".git", "")
    else:
        return repo_url


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
    # cd_repo might not point to the root of the
    # the git repository because git log --find-object works
    # from anywhere in the repo, hence, let's get the root
    git_root = subprocess.run(
        "git rev-parse --show-toplevel",
        shell=True,
        capture_output=True,
        cwd=cd_repo,
    ).stdout.decode()
    command = f"git ls-tree -r {commit_hash} | grep -E {blob_hash}"
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        cwd=git_root,
    )
    print("actual call")
    print(result.stdout.decode())
    if result.returncode != 0 and result.stderr.decode() != "":
        raise RuntimeError(f"{command}\n{result.stderr.decode()}")
    if len(result.stdout.decode()) == 0:
        raise RuntimeError(
            f"Could not find filepath within git repo running:\n{command}"
            f"\nin repo: {git_root}"
        )
    filepath = result.stdout.decode().split()[-1]
    return filepath


def get_transform_reference_from_git_repo(path: Path):
    blob_hash = hash_code(path).hexdigest()
    cd_repo = None
    print(Path.cwd())
    result = get_git_commit_hash(blob_hash, cd_repo=cd_repo)
    commit_hash = result.stdout.decode()
    print(Path.cwd())
    print("with HEAD")
    result = subprocess.run(
        "git ls-tree -r HEAD",
        shell=True,
        capture_output=True,
        cwd=cd_repo,
    )
    print(result.stdout.decode())
    print("with commit hash")
    result = subprocess.run(
        f"git ls-tree -r {commit_hash}",
        shell=True,
        capture_output=True,
        cwd=cd_repo,
    )
    print(result.stdout.decode())
    print(commit_hash, cd_repo, result.returncode)
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
