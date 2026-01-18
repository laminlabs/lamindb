import lamindb as ln
from pyinstrument import Session

threshold = 3.2
session = Session.load("profiling_session.pyisession")
duration = session.duration

ln.connect("laminlabs/lamindata")
ln.track("eraGM939WmQO")
sheet = ln.Record.get(name="import_lamindb.py")
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
