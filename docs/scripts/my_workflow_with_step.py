import lamindb as ln


@ln.step()
def subset_dataframe(
    artifact: ln.Artifact,
    subset_rows: int = 2,
    subset_cols: int = 2,
) -> ln.Artifact:
    df = artifact.load()
    new_data = df.iloc[:subset_rows, :subset_cols]
    new_key = artifact.key.replace(".parquet", "_subsetted.parquet")
    return ln.Artifact.from_dataframe(new_data, key=new_key).save()


@ln.flow()
def ingest_dataset(key: str, subset: bool = False) -> ln.Artifact:
    df = ln.examples.datasets.mini_immuno.get_dataset1()
    artifact = ln.Artifact.from_dataframe(df, key=key).save()
    if subset:
        artifact = subset_dataframe(artifact)
    return artifact


if __name__ == "__main__":
    ingest_dataset(key="my_analysis/dataset.parquet", subset=True)
