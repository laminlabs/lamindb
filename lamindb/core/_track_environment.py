from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
from lamin_utils import logger

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


def _track_pip_freeze(env_dir: Path) -> bool:
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
        )
    except OSError as e:
        logger.warning(f"could not run pip freeze: {e}")
        return False
    if result.returncode != 0:
        err = (result.stderr or "").strip()
        logger.warning(
            "could not run pip freeze"
            + (f": {err}" if err else f" (exit code {result.returncode})")
        )
        return False
    if not result.stdout.strip():
        logger.warning("pip freeze returned empty output, skipping environment tracking")
        return False
    env_dir.mkdir(parents=True, exist_ok=True)
    filepath = env_dir / "run_env_pip.txt"
    filepath.write_text(result.stdout)
    logger.info(f"tracked pip freeze > {filepath}")
    return True


def track_python_environment(run: Run) -> None:
    env_dir = ln_setup.settings.cache_dir / "environments" / f"run_{run.uid}"
    # Prefer pixi lockfile when running inside a pixi environment.
    if _track_pixi_lock(env_dir):
        return
    # we don't create a conda environment.yml mostly for its slowness
    if _track_pip_freeze(env_dir):
        return
    logger.warning(
        "could not track Python environment"
        " (no pixi.lock found, pip freeze unavailable)"
    )
