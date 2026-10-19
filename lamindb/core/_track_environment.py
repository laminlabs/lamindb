from __future__ import annotations

import os
import shutil
import subprocess
import sys
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamindb_setup import logger

if TYPE_CHECKING:
    from lamindb.models import Run


def _find_pixi_project_root() -> Path | None:
    """Return the pixi project root if this process is running inside a pixi env."""
    # sys.prefix inside a pixi env looks like: .../<project>/.pixi/envs/<name>
    prefix = Path(sys.prefix).resolve()
    parts = prefix.parts
    try:
        idx = parts.index(".pixi")
    except ValueError:
        idx = -1
    if idx > 0 and idx + 1 < len(parts) and parts[idx + 1] == "envs":
        return Path(*parts[:idx])
    # Fallback: walk up from cwd looking for pixi.toml
    for path in (Path.cwd(), *Path.cwd().parents):
        if (path / "pixi.toml").exists():
            return path
    return None


def _track_pixi_lock(env_dir: Path) -> bool:
    root = _find_pixi_project_root()
    if root is None:
        return False
    lock = root / "pixi.lock"
    if not lock.is_file() or lock.stat().st_size == 0:
        return False
    env_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(lock, env_dir / "pixi.lock")
    logger.info(f"tracked pixi.lock > {env_dir / 'pixi.lock'}")
    return True


def _pyvenv_created_by_uv(prefix: Path) -> bool:
    """uv writes `uv = <version>` into pyvenv.cfg; the stdlib venv does not."""
    try:
        text = (prefix / "pyvenv.cfg").read_text()
    except (OSError, UnicodeError):
        return False
    for line in text.splitlines():
        key, sep, _ = line.partition("=")
        if sep and key.strip() == "uv":
            return True
    return False


def _uv_lock_candidates() -> list[Path]:
    """Paths that may hold the uv lockfile for this process.

    `uv run` exports `UV` as the path of the uv executable, not the project.
    `Path(UV).parents[2] / "uv.lock"` is the project lock only when that
    executable lives inside the project (`<project>/.venv/bin/uv`). The usual
    layout keeps the lock next to the project environment: `sys.prefix` is
    `<project>/.venv`, so the lock is its parent, which is also
    `Path(sys.executable).parents[2]`.
    """
    candidates: list[Path] = []
    prefix = Path(sys.prefix)
    if prefix.name == ".venv":
        candidates.append(prefix.parent / "uv.lock")
    executable = Path(sys.executable)
    if len(executable.parents) > 2:
        candidates.append(executable.parents[2] / "uv.lock")
    uv = os.environ.get("UV")
    if uv:
        uv_path = Path(uv)
        if len(uv_path.parents) > 2:
            candidates.append(uv_path.parents[2] / "uv.lock")
    # Same fallback as pixi: a lockfile above the working directory.
    candidates.extend(path / "uv.lock" for path in (Path.cwd(), *Path.cwd().parents))
    return candidates


def _find_uv_lock() -> Path | None:
    """Return `uv.lock` if this process is running inside a uv environment."""
    # `uv run` sets UV. Activating a uv virtualenv does not, but uv stamps
    # pyvenv.cfg, and those envs have no pip either.
    if "UV" not in os.environ and not _pyvenv_created_by_uv(Path(sys.prefix)):
        return None
    seen: set[Path] = set()
    for lock in _uv_lock_candidates():
        if lock in seen:
            continue
        seen.add(lock)
        if lock.is_file() and lock.stat().st_size > 0:
            return lock
    return None


def _track_uv_lock(env_dir: Path) -> bool:
    lock = _find_uv_lock()
    if lock is None:
        return False
    env_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(lock, env_dir / "uv.lock")
    logger.info(f"tracked uv.lock > {env_dir / 'uv.lock'}")
    return True


def _track_pip_freeze(env_dir: Path) -> str | None:
    """Snapshot `pip freeze`. Return an error message on failure."""
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
        )
    except OSError as e:
        return f"could not run pip freeze: {e}"
    if result.returncode != 0:
        err = (result.stderr or "").strip()
        detail = f": {err}" if err else f" (exit code {result.returncode})"
        return f"could not run pip freeze{detail}"
    if not result.stdout.strip():
        return "pip freeze returned empty output"
    env_dir.mkdir(parents=True, exist_ok=True)
    filepath = env_dir / "run_env_pip.txt"
    filepath.write_text(result.stdout)
    logger.info(f"tracked pip freeze > {filepath}")
    return None


def track_python_environment(run: Run) -> None:
    env_dir = ln_setup.settings.cache_dir / "environments" / f"run_{run.uid}"
    # Prefer a lockfile when running inside pixi or uv. uv virtualenvs ship
    # without pip, so `python -m pip freeze` fails there.
    if _track_pixi_lock(env_dir) or _track_uv_lock(env_dir):
        return
    # we don't create a conda environment.yml mostly for its slowness
    error = _track_pip_freeze(env_dir)
    if error is None:
        return
    warnings.warn(
        f"could not track Python environment (no pixi.lock or uv.lock found): {error}",
        UserWarning,
        stacklevel=2,
    )
