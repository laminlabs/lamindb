from pathlib import Path

import lamindb as ln

sample_id = ln.Feature(name="sample_id", dtype=str).save()
tissue = ln.Feature(name="tissue", dtype=str).save()
external_schema = ln.Schema(features=[sample_id, tissue]).save()
schema = ln.Schema(
    name="SIToBeCurated",
    slots={"__external__": external_schema},
).save()

path = Path("reads.fastq.gz")
path.write_bytes(b"@read1\nACGT\n+\n!!!!\n")
artifact = ln.Artifact(
    path,
    key="seq/reads.fastq.gz",
    features={"sample_id": "S1", "tissue": "liver"},
    schema=schema,
).save()
artifact.describe()
