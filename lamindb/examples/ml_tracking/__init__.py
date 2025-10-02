"""Examples and utilities for ML tracking frameworks.

.. autosummary::
   :toctree: .

   create_mlflow_schema
   create_wandb_schema
"""

import lamindb as ln
from lamindb.models import Schema

# The following comments apply to all of the following schemas.
# A different design could leverage records instead where ID and name are bundled.
# This could potentially show up nicer in the UI and allow for markdown comments on MLflow experiments or runs.


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


def create_wandb_schema() -> Schema:
    """Creates a W&B schema that covers all W&B project and run related features.

    Saves the following features:
    - wandb_run_id
    - wandb_run_name
    - wandb_project
    - wandb_entity
    - wandb_user
    - wandb_state
    - wandb_url
    - wandb_tags
    - wandb_group
    - wandb_job_type
    - start_time
    - end_time
    """
    wandb_features = [
        ln.Feature(name="wandb_run_id", dtype=str).save(),
        ln.Feature(name="wandb_run_name", dtype=str).save(),
        ln.Feature(name="wandb_project", dtype=str).save(),
        ln.Feature(name="wandb_entity", dtype=str).save(),
        ln.Feature(name="wandb_user", dtype=str).save(),
        ln.Feature(name="wandb_state", dtype=str).save(),
        ln.Feature(name="wandb_url", dtype=str).save(),
        ln.Feature(name="wandb_tags", dtype=str).save(),
        ln.Feature(name="wandb_group", dtype=str).save(),
        ln.Feature(name="wandb_job_type", dtype=str).save(),
        ln.Feature(name="start_time", dtype=int).save(),
        ln.Feature(name="end_time", dtype=int).save(),
    ]
    schema = ln.Schema(
        name="W&B schema",
        features=wandb_features,
        otype="DataFrame",
        minimal_set=False,
        maximal_set=True,
        coerce_dtype=True,
    ).save()
    return schema
