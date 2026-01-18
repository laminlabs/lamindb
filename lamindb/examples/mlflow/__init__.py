"""Examples and utilities for Mlflow.

.. autofunction:: save_mlflow_features
"""

import lamindb as ln


def save_mlflow_features():
    """Saves all MLflow experiment and run related features.

    Saves the following features:

    - mlflow_run_id
    - mlflow_run_name
    - mlflow_experiment_id
    - mlflow_experiment_name
    - mlflow_user_id
    - mlflow_status
    - mlflow_lifecycle_stage
    - mlflow_artifact_uri
    - mlflow_start_time
    - mlflow_end_time
    """
    mlflow_type = ln.Feature(name="MLflow", is_type=True).save()
    ln.Feature(name="mlflow_run_id", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_run_name", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_experiment_id", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_experiment_name", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_user_id", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_status", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_lifecycle_stage", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_artifact_uri", dtype=str, type=mlflow_type).save()
    ln.Feature(name="mlflow_start_time", dtype=int, type=mlflow_type).save()
    ln.Feature(name="mlflow_end_time", dtype=int, type=mlflow_type).save()
