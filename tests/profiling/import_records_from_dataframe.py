from datetime import datetime
from random import Random

import pandas as pd
from laminprofiler import laminprofiler


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


@laminprofiler("EQ44qsBEYgynzl4V")
def profile_import_records_from_dataframe():
    """Time Record.from_dataframe().save() for 100 generated rows on lamindata."""
    import lamindb as ln

    ln.connect("laminlabs/lamindata")
    feature_names = [
        "age_or_mean_of_age_range",
        "array_col",
        "cell_type_by_model",
    ]
    n_rows = 100
    rng = Random(0)
    features = ln.Feature.filter(name__in=feature_names)
    dtypes_by_feature = {feature.name: feature.dtype_as_str for feature in features}
    data: dict[str, list] = {}
    for feature in features:
        data[feature.name] = generate_values(
            dtypes_by_feature[feature.name], n_rows, rng
        )
    df = pd.DataFrame(data)
    type_name = (
        "test-import-records-from-dataframe-"
        f"{datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f')}"
    )
    ln.Record.from_dataframe(df, type=type_name).save()
