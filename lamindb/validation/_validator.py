from typing import Dict, Iterable, Optional

import pandas as pd
from lamin_utils import colors, logger
from lnschema_core.types import FieldAttr

import lamindb as ln

from ._lookup import Lookup
from ._register import register_artifact, register_labels
from ._validate import validate_categories_in_df


class ValidationError(ValueError):
    """Validation error."""

    pass


class Validator:
    """Lamin validator.

    Args:
        df: The DataFrame object to validate.
        fields: A dictionary mapping column to registry_field.
            For example:
            {"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name}
        using: The reference instance containing registries to validate against.
        verbosity: The verbosity level.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        fields: Dict[str, FieldAttr] = None,
        feature_field: FieldAttr = ln.Feature.name,
        using: str = None,
        verbosity: str = "hint",
        **kwargs,
    ) -> None:
        """Validate an AnnData object."""
        if fields is None:
            fields = {}
        self._df = df
        self._fields = fields
        self._feature_field = feature_field
        self._using = using
        ln.settings.verbosity = verbosity
        self._artifact = None
        self._collection = None
        self._validated = False
        self._kwargs: Dict = {}
        self._add_kwargs(**kwargs)
        self.register_features()

    @property
    def fields(self) -> Dict:
        """Return the columns fields to validate against."""
        return self._fields

    def _add_kwargs(self, **kwargs):
        for k, v in kwargs.items():
            self._kwargs[k] = v

    def _register_labels_all(self, validated_only: bool = True, **kwargs):
        """Register labels for all features."""
        for name in self.fields.keys():
            logger.info(f"registering labels for '{name}'")
            self.register_labels(feature=name, validated_only=validated_only, **kwargs)

    def lookup(self, using: Optional[str] = None) -> Lookup:
        """Lookup features and labels.

        Args:
            using: The instance where the lookup is performed.
                if None (default), the lookup is performed on the instance specified in "using" parameter of the Validator.
                if "public", the lookup is performed on the public reference.
        """
        fields = {**{"feature": self._feature_field}, **self.fields}
        return Lookup(fields=fields, using=using or self._using)

    def register_features(self, validated_only: bool = True) -> None:
        """Register features records."""
        # always register features specified as the fields keys
        missing_columns = [i for i in self.fields.keys() if i not in self._df]
        if len(missing_columns) > 0:
            raise ValueError(
                f"columns {missing_columns} are not found in the data object!"
            )
        register_labels(
            values=list(self.fields.keys()),
            field=self._feature_field,
            feature_name="feature",
            using=self._using,
            validated_only=False,
            df=self._df,
            kwargs=self._kwargs,
        )
        # register the rest of the columns based on validated_only
        additional_columns = [
            i for i in self._df.columns if i not in self.fields.keys()
        ]
        if len(additional_columns) > 0:
            register_labels(
                values=additional_columns,
                field=self._feature_field,
                feature_name="feature",
                using=self._using,
                validated_only=validated_only,
                df=self._df,
                kwargs=self._kwargs,
            )

    def register_labels(self, feature: str, validated_only: bool = True, **kwargs):
        """Register labels for a feature.

        Args:
            feature: The name of the feature to register.
            validated_only: Whether to register only validated labels.
            **kwargs: Additional keyword arguments.
        """
        if feature == "all":
            self._register_labels_all(validated_only=validated_only, **kwargs)
        elif feature == "feature":
            self.register_features(validated_only=validated_only)
        else:
            if feature not in self.fields:
                raise ValueError(f"feature {feature} is not part of the fields!")
            register_labels(
                values=self._df[feature].unique().tolist(),
                field=self.fields.get(feature),
                feature_name=feature,
                using=self._using,
                validated_only=validated_only,
                kwargs=kwargs,
            )

    def validate(
        self,
        **kwargs,
    ) -> bool:
        """Validate variables and categorical observations.

        Returns:
            whether the AnnData object is validated
        """
        self._add_kwargs(**kwargs)
        self._validated = validate_categories_in_df(
            self._df,
            fields=self.fields,
            using=self._using,
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
            raise ValidationError(
                f"data object is not validated, please run {colors.yellow('validate()')}!"
            )

        # make sure all labels are registered in the current instance
        verbosity = ln.settings.verbosity
        try:
            ln.settings.verbosity = "warning"
            self.register_labels("all")

            self._artifact = register_artifact(
                self._df,
                description=description,
                fields=self.fields,
                feature_field=self._feature_field,
                **self._kwargs,
            )
        finally:
            ln.settings.verbosity = verbosity

        return self._artifact

    def register_collection(
        self,
        artifact: ln.Artifact | Iterable[ln.Artifact],
        name: str,
        description: Optional[str] = None,
        reference: Optional[str] = None,
        reference_type: Optional[str] = None,
    ) -> ln.Collection:
        """Register a collection from artifact/artifacts.

        Args:
            artifact: one or several registered Artifacts
            name: title of the publication
            description: description of the publication
            reference: accession number (e.g. GSE#, E-MTAB#, etc.)
            reference_type: source type (e.g. GEO, ArrayExpress, SRA, etc.)
        """
        collection = ln.Collection(
            artifact,
            name=name,
            description=description,
            reference=reference,
            reference_type=reference_type,
        )
        slug = ln.setup.settings.instance.slug
        if collection._state.adding:
            collection.save()
            logger.success(f"registered collection in {colors.italic(slug)}")
        else:
            collection.save()
            logger.warning(f"collection already exists in {colors.italic(slug)}!")
        if ln.setup.settings.instance.is_remote:
            logger.print(f"🔗 https://lamin.ai/{slug}/collection/{collection.uid}")
        self._collection = collection
        return collection

    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't register any outputs."""
        if ln.run_context.transform is not None:
            ln.Run.filter(
                transform=ln.run_context.transform, output_artifacts=None
            ).exclude(uid=ln.run_context.run.uid).delete()
