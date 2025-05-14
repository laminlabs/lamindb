from __future__ import annotations

import signal
from typing import Callable


def chain_signal_handler(
    signum: int,
    handler: Callable,
):
    previous_handler = signal.getsignal(signum)

    def composite(signo=None, frame=None):
        previous_handler(signo, frame)
        handler(signo, frame)

    signal.signal(signum, composite)
