import argparse
import lamindb as ln

ln.Param(name="run_workflow_subset", dtype=bool).save()


@ln.tracked()
def subset_dataframe(
    artifact: ln.Artifact,
    subset_rows: int = 2,
    subset_cols: int = 2,
    run: ln.Run | None = None,
) -> ln.Artifact:
    dataset = artifact.load(is_run_input=run)
    new_data = dataset.iloc[:subset_rows, :subset_cols]
    new_key = artifact.key.replace(".parquet", "_subsetted.parquet")
    return ln.Artifact.from_df(new_data, key=new_key, run=run).save()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--subset", action="store_true")
    args = p.parse_args()

    params = {"run_workflow_subset": args.subset}

    ln.track(params=params)

    if args.subset:
        df = ln.core.datasets.small_dataset1(otype="DataFrame")
        artifact = ln.Artifact.from_df(df, key="my_analysis/dataset.parquet").save()
        subsetted_artifact = subset_dataframe(artifact)

    ln.finish()
