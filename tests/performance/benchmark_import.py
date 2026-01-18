from pyinstrument import Profiler

profiler = Profiler()
profiler.start()


profiler.stop()

duration = profiler.session.duration
threshold = 2.5

if duration > 2.5:
    print(f"ERROR: Import time {duration:.3f}s exceeds threshold {threshold:.3f}s")
    raise SystemExit(1)
