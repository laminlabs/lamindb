"""Examples and utilities for MLflow runs.

.. autosummary::
   :toctree: .

   create_mlflow_schema
"""

import lamindb as ln
from lamindb.models import Schema


def create_mlflow_schema() -> Schema:
    """Creates a MLflow schema that covers all MLflow experiment and run related features.

    Saves the following features:

    - mlflow_run_id
    - mlflow_run_name
    - mlflow_experiment_id
    - mlflow_experiment_name
    - mlflow_user_id
    - mlflow_status
    - mlflow_lifecycle_stage
    - mlflow_artifact_uri
    - start_time
    - end_time
    """
    # A different design could leverage records instead where ID and name are bundled.
    # This could potentially show up nicer in the UI and allow for markdown comments on MLflow experiments or runs.
    mlflow_features = [
        ln.Feature(name="mlflow_run_id", dtype=str).save(),
        ln.Feature(name="mlflow_run_name", dtype=str).save(),
        ln.Feature(name="mlflow_experiment_id", dtype=str).save(),
        ln.Feature(name="mlflow_experiment_name", dtype=str).save(),
        ln.Feature(name="mlflow_user_id", dtype=str).save(),
        ln.Feature(name="mlflow_status", dtype=str).save(),
        ln.Feature(name="mlflow_lifecycle_stage", dtype=str).save(),
        ln.Feature(name="mlflow_artifact_uri", dtype=str).save(),
        ln.Feature(name="start_time", dtype=int).save(),
        ln.Feature(name="end_time", dtype=int).save(),
    ]

    schema = ln.Schema(
        name="MLflow schema",
        features=mlflow_features,
        otype="DataFrame",
        minimal_set=False,
        maximal_set=True,
        coerce_dtype=True,
    ).save()

    return schema
