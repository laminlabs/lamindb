import re

import lamindb as ln

threshold = 3.4

# Parse duration from pyinstrument text output
with open("profile.txt") as f:
    content = f.read()
    # Extract duration from line like: "Duration: 2.315     CPU time: 2.216"
    match = re.search(r"Duration:\s+([\d.]+)", content)
    duration = float(match.group(1)) if match else 1.0

print(content)

print(f"Extracted duration: {duration:.3f}s")

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
