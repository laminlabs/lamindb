"""Full/meta-package module for the `lamindb` distribution."""

from __future__ import annotations

import re
from pathlib import Path

_INIT_FILE = Path(__file__).parent / "lamindb" / "__init__.py"
_MATCH = re.search(r'__version__\s*=\s*"([^"]+)"', _INIT_FILE.read_text())
if _MATCH is None:
    raise RuntimeError(f"Could not parse __version__ from {_INIT_FILE}")

__version__ = _MATCH.group(1)
