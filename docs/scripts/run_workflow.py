import lamindb as ln


@ln.flow()
def ingest_dataset(key: str) -> ln.Artifact:
    """Ingest a dataset.

    Args:
        key: The key for the dataset, e.g., "my_analysis/dataset.parquet".
    """
    df = ln.examples.datasets.mini_immuno.get_dataset1()
    artifact = ln.Artifact.from_dataframe(df, key=key).save()
    return artifact


if __name__ == "__main__":
    ingest_dataset(key="my_analysis/dataset.parquet")
