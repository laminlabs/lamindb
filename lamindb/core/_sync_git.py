from __future__ import annotations

import subprocess
from pathlib import Path

from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core.hashing import hash_code

from ..core._settings import sanitize_git_repo_url, settings
from ..errors import BlobHashNotFound


def get_git_repo_from_remote(url: str | None = None, depth: int | None = 10) -> Path:
    """Clone the git repository if not already cloned.

    If `depth` is provided, a shallow clone is performed and no tags are fetched.
    """
    repo_url = url or settings.sync_git_repo
    repo_dir = setup_settings.cache_dir / repo_url.split("/")[-1]
    if repo_dir.exists():
        logger.debug(f"git repo {repo_dir} already exists locally")
        return repo_dir
    logger.important(
        f"running outside of synched git repo, cloning {repo_url} into {repo_dir}"
    )
    args = ["git", "clone", f"{repo_url}.git"]
    if depth is not None:
        # if depth is provided, will not fetch tags
        args += ["--depth", f"{depth}"]
    result = subprocess.run(
        args,
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
            f"❌ Did not find blob hash {blob_hash} in git repo: {settings.sync_git_repo}\n"
            f"Did you commit & push the script to the remote repo? -> {path}"
        )
    gitpath = get_filepath_within_git_repo(commit_hash, blob_hash, repo_dir)
    reference = f"{settings.sync_git_repo}/blob/{commit_hash}/{gitpath}"
    return reference


def get_and_validate_git_metadata(
    url: str,
    path: str,
    version: str | None = None,
    branch: str | None = None,
) -> tuple[str, str]:
    """Get metadata from a git repository.

    Args:
        url: Git repository URL (e.g., "https://github.com/user/repo")
        path: Path to the main script within the repository
        version: Optional version/tag to checkout
        branch: Optional branch name (defaults to repository's default branch)

    Returns:
        Dictionary containing:
            - commit_hash: The current commit hash
            - url: The repository URL
            - main_script: Path to the main script
            - revision: The version/tag (if provided)
            - branch: The branch name

    Raises:
        RuntimeError: If git operations fail
        FileNotFoundError: If the specified path does not exist in the repository
    """
    url = sanitize_git_repo_url(url)
    repo_dir = get_git_repo_from_remote(url, depth=None)

    # Determine the branch to use
    if branch is None:
        # Get the default branch if not specified
        result_str = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "origin/HEAD"],
            capture_output=True,
            cwd=repo_dir,
            text=True,
        )
        if result_str.returncode == 0:
            branch = result_str.stdout.strip().split("/")[-1]
        else:
            branch = "main"  # fallback to main

    # Fetch the latest changes
    subprocess.run(
        ["git", "fetch", "origin"],
        capture_output=True,
        cwd=repo_dir,
        check=True,
    )

    # Checkout the specified version or branch
    if version is not None:
        # Version takes precedence - checkout the tag/version
        result = subprocess.run(
            ["git", "checkout", version],
            capture_output=True,
            cwd=repo_dir,
        )
        if result.returncode != 0:
            raise ValueError(
                f"Failed to checkout version {version}: {result.stderr.decode()}"
            )
        logger.info(f"checked out version {version}")
    else:
        # Checkout the branch
        result = subprocess.run(
            ["git", "checkout", f"origin/{branch}"],
            capture_output=True,
            cwd=repo_dir,
        )
        if result.returncode != 0:
            raise ValueError(
                f"Failed to checkout branch {branch}: {result.stderr.decode()}"
            )
        logger.info(f"checked out branch {branch}")

    # Get the current commit hash
    result_str = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True,
        cwd=repo_dir,
        text=True,
    )
    if result_str.returncode != 0:
        raise RuntimeError(f"Failed to get commit hash: {result_str.stderr}")

    commit_hash = result_str.stdout.strip()

    assert (  # noqa: S101
        len(commit_hash) == 40
    ), f"commit hash |{commit_hash}| is not 40 characters long"

    # Verify that the path exists as a file in the repository
    file_path = repo_dir / path
    if not file_path.exists():
        raise FileNotFoundError(f"Path '{path}' does not exist in repository {url}")
    if not file_path.is_file():
        raise FileNotFoundError(
            f"Path '{path}' exists but is not a file in repository {url}"
        )
    return url, commit_hash
