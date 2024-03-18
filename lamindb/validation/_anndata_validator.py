from typing import Dict, Optional

import anndata as ad
from lnschema_core.types import FieldAttr
from pandas.core.api import DataFrame as DataFrame

import lamindb as ln

from ._lookup import Lookup
from ._register import register_artifact, register_labels
from ._validate import validate_anndata
from ._validator import Validator


class AnnDataValidator(Validator):
    """Lamin AnnData validator.

    Args:
        adata: The AnnData object to validate.
        var_field: The registry field to validate variables index against.
        obs_fields: A dictionary mapping obs_column to registry_field.
            For example:
            {"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name}
        using: The reference instance containing registries to validate against.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        var_field: FieldAttr,
        obs_fields: Dict[str, FieldAttr],
        using: str = "default",
        verbosity: str = "hint",
        **kwargs,
    ) -> None:
        self._adata = adata
        super().__init__(
            df=self._adata.obs,
            fields=obs_fields,
            using=using,
            verbosity=verbosity,
            **kwargs,
        )
        self._obs_fields = obs_fields
        self._var_field = var_field
        self._fields = {"variables": var_field, **obs_fields}

    @property
    def var_field(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def obs_fields(self) -> Dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(self, using: Optional[str] = None) -> Lookup:
        """Lookup features and labels."""
        fields = {
            **{"feature": ln.Feature.name, "variables": self.var_field},
            **self.obs_fields,
        }
        return Lookup(fields=fields, using=using or self._using)

    def _register_variables(self, validated_only: bool = True, **kwargs):
        """Register variable records."""
        self._add_kwargs(**kwargs)
        register_labels(
            values=self._adata.var_names,
            field=self.var_field,
            feature_name="variables",
            using=self._using,
            validated_only=validated_only,
            kwargs=self._kwargs,
        )

    def register_labels(self, feature: str, validated_only: bool = True, **kwargs):
        """Register labels for the given feature.

        Args:
            feature: The feature to register labels for.
                if "variables", register variables.
            validated_only: If True, only register validated labels.
            **kwargs: Additional metadata needed.
        """
        if feature == "variables":
            self._register_variables(validated_only=validated_only, **kwargs)
        else:
            super().register_labels(feature, validated_only, **kwargs)

    def validate(self, **kwargs) -> bool:
        """Validate variables and categorical observations."""
        self._add_kwargs(**kwargs)
        self._validated = validate_anndata(
            self._adata,
            var_field=self.var_field,
            obs_fields=self.obs_fields,
            **self._kwargs,
        )

        return self._validated

    def register_artifact(
        self,
        description: str,
        **kwargs,
    ) -> ln.Artifact:
        """Register the validated AnnData and metadata.

        Args:
            description: description of the AnnData object
            **kwargs: object level metadata

        Returns:
            a registered artifact record
        """
        self._add_kwargs(**kwargs)
        if not self._validated:
            raise ValueError("please run `validate()` first!")

        self._artifact = register_artifact(
            self._adata,
            description=description,
            var_field=self.var_field,
            fields=self.obs_fields,
            **self._kwargs,
        )

        return self._artifact
