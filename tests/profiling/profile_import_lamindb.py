from pyinstrument import Profiler

profiler = Profiler()
profiler.start()

import lamindb as ln  # noqa: E402

profiler.stop()

duration = profiler.last_session.duration
threshold = 3.2

print(profiler.output_text())

ln.connect("laminlabs/lamindata")
ln.track()
sheet = ln.Record.get(name="ImportLaminDB")
record = ln.Record(type=sheet).save()
record.features.add_values(
    {
        "duration_in_sec": duration,
        "lamindb_version": ln.__version__,
    }
)
ln.finish()

if duration > threshold:
    print(f"ERROR: Import time {duration:.3f}s exceeds threshold {threshold:.3f}s")
    raise SystemExit(1)
