"""Example schemas.

.. autofunction:: valid_features
.. autofunction:: anndata_ensembl_gene_ids_and_valid_features_in_obs
.. autofunction:: a_feature_type
.. autofunction:: anndata_ensembl_gene_ids_and_obs_feature_type

"""

from ._anndata import (
    anndata_ensembl_gene_ids_and_obs_feature_type,
    anndata_ensembl_gene_ids_and_valid_features_in_obs,
)
from ._simple import a_feature_type, valid_features
