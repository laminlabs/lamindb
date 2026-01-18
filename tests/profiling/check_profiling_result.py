import lamindb as ln

threshold = 3.2
duration = 1

ln.connect("laminlabs/lamindata")
ln.track("eraGM939WmQO")
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
