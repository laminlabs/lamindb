from __future__ import annotations

import subprocess
from pathlib import Path

from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.hashing import hash_code

from ..core._settings import sanitize_git_repo_url, settings


class BlobHashNotFound(SystemExit):
    pass


def get_git_repo_from_remote() -> Path:
    repo_url = settings.sync_git_repo
    repo_dir = setup_settings.cache_dir / repo_url.split("/")[-1]
    if repo_dir.exists():
        logger.warning(f"git repo {repo_dir} already exists locally")
        return repo_dir
    logger.important(
        f"running outside of synched git repo, cloning {repo_url} into {repo_dir}"
    )
    result = subprocess.run(
        ["git", "clone", "--depth", "10", f"{repo_url}.git"],
        capture_output=True,
        cwd=setup_settings.cache_dir,
    )
    if result.returncode != 0 or not repo_dir.exists():
        raise RuntimeError(result.stderr.decode())
    return repo_dir


def check_local_git_repo() -> bool:
    result = subprocess.run(
        ["git", "config", "--get remote.origin.url"],
        capture_output=True,
    )
    result_str = result.stdout.decode().strip()
    if result_str == "":
        # running-not-in-a-git-repo
        return False
    else:
        remote_url = sanitize_git_repo_url(result_str)
        if remote_url == settings.sync_git_repo:
            # running-in-correct-git-repo
            return True
        else:
            # running-outside-of-correct-git-repo
            return False


def get_git_commit_hash(blob_hash: str, repo_dir: Path | None = None) -> str | None:
    # Fetch all remote branches so that we can also search them
    fetch_command = ["git", "fetch", "origin", "+refs/heads/*:refs/remotes/origin/*"]
    subprocess.run(fetch_command, cwd=repo_dir, check=True)

    # Find the commit containing the blob hash in all branches
    command = [
        "git",
        "log",
        "--all",
        f"--find-object={blob_hash}",
        "--pretty=format:%H",
    ]
    result = subprocess.run(
        command,
        capture_output=True,
        cwd=repo_dir,
    )
    # We just care to find one commit
    # Hence, we split by new line ("\n") and use the first one
    commit_hash = result.stdout.decode().split("\n")[0]

    if not commit_hash or result.returncode == 1:
        return None

    default_branch = (
        subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "origin/HEAD"],
            capture_output=True,
            cwd=repo_dir,
            text=True,
        )
        .stdout.strip()
        .split("/")[-1]
    )

    # Find all branches containing the commit
    commit_containing_branches = subprocess.run(
        ["git", "branch", "--all", "--contains", commit_hash],
        capture_output=True,
        cwd=repo_dir,
        text=True,
    ).stdout.split("\n")

    # Clean up branch names and filter out the default branch
    commit_containing_branches = [
        branch.strip().replace("remotes/", "")
        for branch in commit_containing_branches
        if branch.strip()
    ]
    non_default_branches = [
        branch for branch in commit_containing_branches if default_branch not in branch
    ]

    if non_default_branches:
        logger.warning(
            f"code blob hash {blob_hash} was found in non-default branch(es): {', '.join(non_default_branches)}"
        )

    assert (  # noqa: S101
        len(commit_hash) == 40
    ), f"commit hash |{commit_hash}| is not 40 characters long"

    return commit_hash


def get_filepath_within_git_repo(
    commit_hash: str, blob_hash: str, repo_dir: Path | None
) -> str:
    # repo_dir might not point to the root of the
    # the git repository because git log --find-object works
    # from anywhere in the repo, hence, let's get the root
    repo_root = (
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            cwd=repo_dir,
        )
        .stdout.decode()
        .strip()
    )
    # Run the git commands separately to circumvent spawning a shell
    git_command = ["git", "ls-tree", "-r", commit_hash]
    git_process = subprocess.Popen(
        git_command,
        stdout=subprocess.PIPE,
        cwd=repo_root,
    )

    grep_command = ["grep", "-E", blob_hash]
    result = subprocess.run(
        grep_command,
        stdin=git_process.stdout,
        capture_output=True,
        cwd=repo_root,
    )

    # Close the stdout to allow git_process to receive a SIGPIPE if grep_command exits
    git_process.stdout.close()
    git_process.wait()

    command = " ".join(git_command) + " | " + " ".join(grep_command)
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
    commit_hash = None
    if check_local_git_repo():
        repo_dir = None
    else:
        repo_dir = get_git_repo_from_remote()
    commit_hash = get_git_commit_hash(blob_hash, repo_dir=repo_dir)
    if commit_hash is None:
        if repo_dir is None:
            repo_dir = Path.cwd()
        raise BlobHashNotFound(
            f"❌ Did not find blob hash {blob_hash} in git repo ({settings.sync_git_repo}) {repo_dir}\n"
            f"Did you commit the script? -> {path}"
        )
    gitpath = get_filepath_within_git_repo(commit_hash, blob_hash, repo_dir)
    reference = f"{settings.sync_git_repo}/blob/{commit_hash}/{gitpath}"
    return reference
