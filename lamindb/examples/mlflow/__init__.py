"""Examples and utilities for MLflow runs.

.. autosummary::
   :toctree: .

   save_default_values

"""

import lamindb as ln


def save_default_values() -> None:
    """Save commonly tracked default values of MLflow runs to the instance.

    Stores the following features:

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
    ln.Feature(name="mlflow_run_id", dtype=str).save()
    ln.Feature(name="mlflow_run_name", dtype=str).save()
    ln.Feature(name="mlflow_experiment_id", dtype=str).save()
    ln.Feature(name="mlflow_experiment_name", dtype=str).save()
    ln.Feature(name="mlflow_user_id", dtype=str).save()
    ln.Feature(name="mlflow_status", dtype=str).save()
    ln.Feature(name="mlflow_lifecycle_stage", dtype=str).save()
    ln.Feature(name="mlflow_artifact_uri", dtype=str).save()
    ln.Feature(name="start_time", dtype=int).save()
    ln.Feature(name="end_time", dtype=int).save()
