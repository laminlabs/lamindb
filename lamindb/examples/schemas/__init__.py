"""Example schemas.

.. autofunction:: valid_features
.. autofunction:: anndata_ensembl_gene_ids_and_valid_features_in_obs

"""

from ._anndata import (
    _anndata_ensembl_gene_ids_and_obs_feature_type,
    anndata_ensembl_gene_ids_and_valid_features_in_obs,
)
from ._simple import _a_feature_type, valid_features
