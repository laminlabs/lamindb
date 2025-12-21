import click
import lamindb as ln


@click.command()
@click.option("--key", required=True)
@ln.flow()
def ingest_dataset(key: str) -> ln.Artifact:
    """Ingest a dataset into LaminDB.

    Args:
        key: The key of the dataset.
    """
    df = ln.examples.datasets.mini_immuno.get_dataset2()
    artifact = ln.Artifact.from_dataframe(df, key=key).save()
    return artifact


if __name__ == "__main__":
    ingest_dataset()
