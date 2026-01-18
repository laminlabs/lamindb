from pyinstrument import Profiler

profiler = Profiler()
profiler.start()

import lamindb  # noqa: E402, F401

profiler.stop()

duration = profiler.last_session.duration
threshold = 2.5

print(profiler.output_text())

if duration > threshold:
    print(f"ERROR: Import time {duration:.3f}s exceeds threshold {threshold:.3f}s")
    raise SystemExit(1)
