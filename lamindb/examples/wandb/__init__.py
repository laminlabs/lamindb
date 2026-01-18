"""Examples and utilities for Weights & Biases.

.. autofunction:: save_wandb_features
"""

import lamindb as ln


def save_wandb_features():
    """Saves all Weights & Biases project and run related features.

    Saves the following features:

    - wandb_run_id
    - wandb_run_name
    - wandb_run_entity
    - wandb_project
    - wandb_state
    - wandb_url
    - wandb_tags
    - wandb_group
    - wandb_job_type
    - timestamp
    - runtime
    """
    wandb_type = ln.Feature(name="Weights & Biases", is_type=True).save()
    ln.Feature(name="wandb_run_id", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_run_name", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_run_entity", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_project", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_state", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_url", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_tags", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_group", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_job_type", dtype=str, type=wandb_type).save()
    ln.Feature(name="wandb_timestamp", dtype=float, type=wandb_type).save()
    ln.Feature(name="wandb_runtime", dtype=float, type=wandb_type).save()
