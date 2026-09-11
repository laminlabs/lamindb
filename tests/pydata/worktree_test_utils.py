from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import lamindb_setup as ln_setup


def run_worktree_case(case: str, tmp_path: Path) -> None:
    settings_home = tmp_path / "settings"
    shutil.copytree(ln_setup.settings.settings_dir, settings_home / ".lamin")
    env = os.environ.copy()
    env["LAMIN_SETTINGS_DIR"] = str(settings_home)
    runner = Path(__file__).with_name("worktree_test_cases.py")
    result = subprocess.run(
        [sys.executable, str(runner), case, str(tmp_path / "worktree-case")],
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, f"stdout: {result.stdout}\nstderr: {result.stderr}"
