import click
import lamindb as ln


@click.command()
@click.option(
    "--key",
    default="my_analysis/dataset.parquet",
    help="The key of the dataset to ingest",
)
@ln.flow()
def ingest_dataset(key: str) -> ln.Artifact:
    """Ingest a dataset into LaminDB.

    Args:
        key: The key of the dataset to ingest, e.g., "my_analysis/dataset.parquet".

    Returns:
        The artifact of the ingested dataset.
    """
    df = ln.examples.datasets.mini_immuno.get_dataset1()
    artifact = ln.Artifact.from_dataframe(df, key=key).save()
    return artifact


if __name__ == "__main__":
    ingest_dataset()
