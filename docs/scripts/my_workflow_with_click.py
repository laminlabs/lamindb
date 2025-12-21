import click
import lamindb as ln


@click.command()
@click.option("--key", required=True)
@ln.flow()
def main(key: str):
    df = ln.examples.datasets.mini_immuno.get_dataset2()
    ln.Artifact.from_dataframe(df, key=key).save()


if __name__ == "__main__":
    main()
