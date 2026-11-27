import argparse
from datetime import datetime
from random import Random
from time import perf_counter

import lamindb as ln
import pandas as pd


def generate_values(dtype: str, n_rows: int, rng: Random):
    cell_types = [
        "T cell",
        "B cell",
        "natural killer cell",
        "monocyte",
        "epithelial cell",
    ]
    if dtype in {"float", "num"}:
        return [round(rng.uniform(0.0, 100.0), 3) for _ in range(n_rows)]
    if dtype.startswith("cat["):
        return [rng.choice(cell_types) for _ in range(n_rows)]
    raise ValueError(f"Unsupported dtype: {dtype}")


@ln.flow("JuJZZEsit1KV")
def main(n_rows: int):
    feature_names = [
        "age_or_mean_of_age_range",
        "array_col",
        "cell_type_by_model",
    ]
    rng = Random(0)
    features = ln.Feature.filter(name__in=feature_names)
    dtypes_by_feature = {feature.name: feature.dtype_as_str for feature in features}

    data: dict[str, list] = {}
    print("Generating random dataframe values...")
    for feature in features:
        data[feature.name] = generate_values(
            dtypes_by_feature[feature.name], n_rows, rng
        )
    df = pd.DataFrame(data)
    print(df.head(5))

    print("Running Record.from_dataframe()...")
    from_dataframe_start = perf_counter()
    records = ln.Record.from_dataframe(
        df,
        type=f"test-import-records-from-dataframe-{datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}",
    )
    from_dataframe_duration_sec = perf_counter() - from_dataframe_start
    print(f"... completed in {from_dataframe_duration_sec:.6f}s")

    print("Saving records...")
    save_start = perf_counter()
    records.save()
    save_duration_sec = perf_counter() - save_start
    print(f"... completed in {save_duration_sec:.6f}s")

    run = ln.context.run
    params = run.params or {}
    params.update(
        {
            "from_dataframe_duration_sec": round(from_dataframe_duration_sec, 6),
            "save_duration_sec": round(save_duration_sec, 6),
        }
    )
    run.params = params
    run.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare and optionally save test Records rows via Record.from_dataframe()."
    )
    parser.add_argument("--rows", type=int, default=100)
    args = parser.parse_args()
    ln.connect("laminlabs/lamindata")
    main(n_rows=args.rows)
