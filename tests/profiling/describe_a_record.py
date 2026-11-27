import lamindb as ln
from laminprofiler import laminprofiler


@laminprofiler("Om6e5GGH77ADM9Yr", total_seconds=20)
def profile_describe_a_record(benchmark):
    """Time loading one lamindata record and calling describe()."""

    def describe():
        db = ln.DB("laminlabs/lamindata")
        record = db.Record.get("mNDJgWFrkWQVW3ox")
        record.describe()

    benchmark.pedantic(describe, rounds=1, iterations=1, warmup_rounds=1)
