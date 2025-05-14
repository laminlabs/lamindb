from __future__ import annotations

import atexit
import signal
import threading
from typing import TYPE_CHECKING

from ..core._signals import chain_signal_handler
from ..core.storage.paths import delete_storage

if TYPE_CHECKING:
    from upath import UPath

cleanup_paths: dict[str, UPath] = {}


def register_cleanup_path(uid: str, path: UPath):
    cleanup_paths[uid] = path


def unregister_cleanup_path(uid: str):
    cleanup_paths.pop(uid, None)


def cleanup(signo=None, frame=None):
    for path in cleanup_paths.values():
        try:
            delete_storage(path, raise_file_not_found_error=False)
        except:  # noqa: S110, E722
            pass


atexit.register(cleanup)

if threading.current_thread() == threading.main_thread():
    chain_signal_handler(signal.SIGTERM, cleanup)
    chain_signal_handler(signal.SIGINT, cleanup)
