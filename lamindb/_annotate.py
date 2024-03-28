from typing import Dict, Iterable, List, Optional, Tuple, Union

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import colors, logger
from lnschema_core import Artifact, Collection, Feature, Registry, Run, ULabel
from lnschema_core.types import FieldAttr


class ValidationError(ValueError):
    """Validation error."""

    pass


class AnnotateLookup:
    """Lookup features and labels from the reference instance."""

    def __init__(
        self, fields: Dict[str, FieldAttr], using: Optional[str] = None
    ) -> None:
        self._fields = fields
        self._using = None if using == "default" else using
        self._using_name = using or ln_setup.settings.instance.slug
        logger.debug(f"Lookup objects from the {colors.italic(self._using_name)}")

    def __getitem__(self, name):
        if name in self._fields:
            registry = self._fields[name].field.model
            if self._using == "public":
                return registry.public().lookup()
            else:
                return get_registry_instance(registry, self._using).lookup()
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def __repr__(self) -> str:
        if len(self._fields) > 0:
            fields = "\n ".join([str([key]) for key in self._fields.keys()])
            return (
                f"Lookup objects from the {colors.italic(self._using_name)}:\n {colors.green(fields)}\n\n"
                "Example:\n    → categories = validator.lookup().['cell_type']\n"
                "    → categories.alveolar_type_1_fibroblast_cell"
            )
        else:
            return colors.warning("No fields are found!")


class DataFrameAnnotator:
    """Annotation flow for a DataFrame object.

    Args:
        df: The DataFrame object to annotate.
        fields: A dictionary mapping column to registry_field.
            For example:
            {"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name}
        feature_field: The field attribute for the feature column.
        using: The reference instance containing registries to validate against.
        verbosity: The verbosity level.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        fields: Optional[Dict[str, FieldAttr]] = None,
        feature_field: FieldAttr = Feature.name,
        using: Optional[str] = None,
        verbosity: str = "hint",
        **kwargs,
    ) -> None:
        from lamindb.core._settings import settings

        self._df = df
        self._fields = fields or {}
        self._feature_field = feature_field
        self._using = using
        settings.verbosity = verbosity
        self._artifact = None
        self._collection = None
        self._validated = False
        self._kwargs: Dict = kwargs
        self.register_features()

    @property
    def fields(self) -> Dict:
        """Return the columns fields to validate against."""
        return self._fields

    def lookup(self, using: Optional[str] = None) -> AnnotateLookup:
        """Lookup features and labels.

        Args:
            using: The instance where the lookup is performed.
                if None (default), the lookup is performed on the instance specified in "using" parameter of the validator.
                if "public", the lookup is performed on the public reference.
        """
        fields = {**{"feature": self._feature_field}, **self.fields}
        return AnnotateLookup(fields=fields, using=using or self._using)

    def register_features(self, validated_only: bool = True) -> None:
        """Register features records."""
        missing_columns = set(self.fields.keys()) - set(self._df.columns)
        if missing_columns:
            raise ValueError(
                f"Columns {missing_columns} are not found in the data object!"
            )

        # Always register features specified as the fields keys
        update_registry(
            values=list(self.fields.keys()),
            field=self._feature_field,
            feature_name="feature",
            using=self._using,
            validated_only=False,
            kwargs=self._kwargs,
        )

        # Register the rest of the columns based on validated_only
        additional_columns = set(self._df.columns) - set(self.fields.keys())
        if additional_columns:
            update_registry(
                values=list(additional_columns),
                field=self._feature_field,
                feature_name="feature",
                using=self._using,
                validated_only=validated_only,
                df=self._df,  # Get the Feature type from df
                kwargs=self._kwargs,
            )

    def update_registry(self, feature: str, validated_only: bool = True, **kwargs):
        """Register labels for a feature.

        Args:
            feature: The name of the feature to register.
            validated_only: Whether to register only validated labels.
            **kwargs: Additional keyword arguments.
        """
        if feature == "all":
            self._update_registry_all(validated_only=validated_only, **kwargs)
        elif feature == "feature":
            self.register_features(validated_only=validated_only)
        else:
            if feature not in self.fields:
                raise ValueError(f"Feature {feature} is not part of the fields!")
            update_registry(
                values=self._df[feature].unique().tolist(),
                field=self.fields[feature],
                feature_name=feature,
                using=self._using,
                validated_only=validated_only,
                kwargs=kwargs,
            )

    def _update_registry_all(self, validated_only: bool = True, **kwargs):
        """Register labels for all features."""
        for name in self.fields.keys():
            logger.info(f"registering labels for '{name}'")
            self.update_registry(feature=name, validated_only=validated_only, **kwargs)

    def validate(self, **kwargs) -> bool:
        """Validate variables and categorical observations.

        Returns:
            Whether the DataFrame is validated.
        """
        self._kwargs.update(kwargs)
        self._validated = validate_categories_in_df(
            self._df,
            fields=self.fields,
            using=self._using,
            **self._kwargs,
        )
        return self._validated

    def register_artifact(self, description: str, **kwargs) -> Artifact:
        """Register the validated DataFrame and metadata.

        Args:
            description: Description of the DataFrame object.
            **kwargs: Object level metadata.

        Returns:
            A registered artifact record.
        """
        from lamindb.core._settings import settings

        self._kwargs.update(kwargs)
        if not self._validated:
            raise ValidationError(
                f"Data object is not validated, please run {colors.yellow('validate()')}!"
            )

        # Make sure all labels are registered in the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            self.update_registry("all")

            self._artifact = register_artifact(
                self._df,
                description=description,
                fields=self.fields,
                feature_field=self._feature_field,
                **self._kwargs,
            )
        finally:
            settings.verbosity = verbosity

        return self._artifact

    def register_collection(
        self,
        artifact: Artifact | Iterable[Artifact],
        name: str,
        description: Optional[str] = None,
        reference: Optional[str] = None,
        reference_type: Optional[str] = None,
    ) -> Collection:
        """Register a collection from artifact/artifacts.

        Args:
            artifact: One or several registered Artifacts.
            name: Title of the publication.
            description: Description of the publication.
            reference: Accession number (e.g. GSE#, E-MTAB#, etc.).
            reference_type: Source type (e.g. GEO, ArrayExpress, SRA, etc.).
        """
        collection = Collection(
            artifact,
            name=name,
            description=description,
            reference=reference,
            reference_type=reference_type,
        )
        slug = ln_setup.settings.instance.slug
        if collection._state.adding:
            collection.save()
            logger.success(f"registered collection in {colors.italic(slug)}")
        else:
            collection.save()
            logger.warning(f"collection already exists in {colors.italic(slug)}!")
        if ln_setup.settings.instance.is_remote:
            logger.print(f"🔗 https://lamin.ai/{slug}/collection/{collection.uid}")
        self._collection = collection
        return collection

    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't register any outputs."""
        from lamindb.core._run_context import run_context

        if run_context.transform is not None:
            Run.filter(transform=run_context.transform, output_artifacts=None).exclude(
                uid=run_context.run.uid
            ).delete()


class AnnDataAnnotator(DataFrameAnnotator):
    """Annotation flow for an AnnData object.

    Args:
        adata: The AnnData object to annotate.
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
        self._var_field = var_field
        super().__init__(
            df=self._adata.obs,
            fields=obs_fields,
            using=using,
            verbosity=verbosity,
            **kwargs,
        )
        self._obs_fields = obs_fields
        self._register_variables()

    @property
    def var_field(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def obs_fields(self) -> Dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(self, using: Optional[str] = None) -> AnnotateLookup:
        """Lookup features and labels."""
        fields = {
            **{"feature": Feature.name, "variables": self.var_field},
            **self.obs_fields,
        }
        return AnnotateLookup(fields=fields, using=using or self._using)

    def _register_variables(self, validated_only: bool = True, **kwargs):
        """Register variable records."""
        self._kwargs.update(kwargs)
        update_registry(
            values=self._adata.var_names,
            field=self.var_field,
            feature_name="variables",
            using=self._using,
            validated_only=validated_only,
            kwargs=self._kwargs,
        )

    def validate(self, **kwargs) -> bool:
        """Validate variables and categorical observations."""
        self._kwargs.update(kwargs)
        self._validated = validate_anndata(
            self._adata,
            var_field=self.var_field,
            obs_fields=self.obs_fields,
            **self._kwargs,
        )
        return self._validated

    def update_registry(self, feature: str, validated_only: bool = True, **kwargs):
        """Register labels for a feature."""
        if feature == "variables":
            self._register_variables(validated_only=validated_only, **kwargs)
        else:
            super().update_registry(feature, validated_only, **kwargs)

    def register_artifact(self, description: str, **kwargs) -> Artifact:
        """Register the validated AnnData and metadata.

        Args:
            description: Description of the AnnData object.
            **kwargs: Object level metadata.

        Returns:
            A registered artifact record.
        """
        self._kwargs.update(kwargs)
        if not self._validated:
            raise ValidationError("Please run `validate()` first!")

        self._artifact = register_artifact(
            self._adata,
            description=description,
            feature_field=self.var_field,
            fields=self.obs_fields,
            **self._kwargs,
        )
        return self._artifact


class Annotate:
    """Annotation flow."""

    @classmethod
    def from_df(
        cls,
        df: pd.DataFrame,
        fields: Optional[Dict[str, FieldAttr]] = None,
        feature_field: FieldAttr = Feature.name,
        using: Optional[str] = None,
        verbosity: str = "hint",
        **kwargs,
    ) -> DataFrameAnnotator:
        return DataFrameAnnotator(
            df=df,
            fields=fields,
            feature_field=feature_field,
            using=using,
            verbosity=verbosity,
            **kwargs,
        )

    @classmethod
    def from_anndata(
        cls,
        adata: ad.AnnData,
        var_field: FieldAttr,
        obs_fields: Dict[str, FieldAttr],
        using: str = "default",
        verbosity: str = "hint",
        **kwargs,
    ) -> AnnDataAnnotator:
        return AnnDataAnnotator(
            adata=adata,
            var_field=var_field,
            obs_fields=obs_fields,
            using=using,
            verbosity=verbosity,
            **kwargs,
        )


def get_registry_instance(registry: Registry, using: Optional[str] = None) -> Registry:
    """Get a registry instance using a specific instance."""
    if using is not None and using != "default":
        return registry.using(using)
    return registry


def standardize_and_inspect(
    values: Iterable[str], field: FieldAttr, registry: Registry, **kwargs
):
    """Standardize and inspect values using a registry."""
    if hasattr(registry, "standardize"):
        values = registry.standardize(values, field=field, mute=True, **kwargs)
    return registry.inspect(values, field=field, mute=True, **kwargs)


def check_registry_organism(
    registry: Registry, organism: Optional[str] = None
) -> Optional[str]:
    """Check if a registry needs an organism and return the organism name."""
    if hasattr(registry, "organism_id"):
        import bionty as bt

        if organism is None and bt.settings.organism is None:
            raise ValueError(
                f"{registry.__name__} registry requires an organism!\n"
                "      → please pass an organism name via organism="
            )
        return organism or bt.settings.organism.name
    return None


def validate_categories(
    values: Iterable[str],
    field: FieldAttr,
    feature_name: str,
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Validate ontology terms in a pandas series using LaminDB registries."""
    from lamindb._from_values import _print_values

    model_field = f"{field.field.model.__name__}.{field.field.name}"
    logger.indent = ""
    logger.info(
        f"inspecting '{colors.bold(feature_name)}' by {colors.italic(model_field)}"
    )
    logger.indent = "   "

    registry = field.field.model
    filter_kwargs = {}
    organism = check_registry_organism(registry, kwargs.get("organism"))
    if organism is not None:
        filter_kwargs["organism"] = organism

    # Inspect the default instance
    inspect_result = standardize_and_inspect(
        values=values, field=field, registry=registry, **filter_kwargs
    )
    non_validated = inspect_result.non_validated

    if using is not None and using != "default" and non_validated:
        registry = get_registry_instance(registry, using)
        # Inspect the using instance
        inspect_result = standardize_and_inspect(
            values=non_validated, field=field, registry=registry, **filter_kwargs
        )
        non_validated = inspect_result.non_validated

    n_non_validated = len(non_validated)
    if n_non_validated == 0:
        logger.success(f"all {feature_name}s are validated")
        return True
    else:
        are = "are" if n_non_validated > 1 else "is"
        print_values = _print_values(non_validated)
        feature_name_print = f".update_registry('{feature_name}')"
        warning_message = (
            f"{colors.yellow(f'{n_non_validated} terms')} {are} not validated: "
            f"{colors.yellow(print_values)}\n      → register terms via "
            f"{colors.yellow(feature_name_print)}"
        )
        logger.warning(warning_message)
        logger.indent = ""
        return False


def validate_categories_in_df(
    df: pd.DataFrame,
    fields: Dict[str, FieldAttr],
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Validate categories in DataFrame columns using LaminDB registries."""
    validated = True
    for feature_name, field in fields.items():
        validated &= validate_categories(
            df[feature_name],
            field=field,
            feature_name=feature_name,
            using=using,
            **kwargs,
        )
    return validated


def validate_anndata(
    adata: ad.AnnData,
    var_field: FieldAttr,
    obs_fields: Dict[str, FieldAttr],
    using: Optional[str] = None,
    **kwargs,
) -> bool:
    """Inspect metadata in an AnnData object using LaminDB registries."""
    if using is not None and using != "default":
        logger.important(
            f"validating metadata using registries of instance {colors.italic(using)}"
        )

    validated_var = validate_categories(
        adata.var.index,
        field=var_field,
        feature_name="variables",
        using=using,
        **kwargs,
    )
    validated_obs = validate_categories_in_df(
        adata.obs, fields=obs_fields, using=using, **kwargs
    )
    return validated_var and validated_obs


def register_artifact(
    data: Union[pd.DataFrame, ad.AnnData],
    description: str,
    fields: Dict[str, FieldAttr],
    feature_field: FieldAttr,
    **kwargs,
) -> Artifact:
    """Register all metadata with an Artifact.

    Args:
        data: The DataFrame or AnnData object to register.
        description: A description of the artifact.
        fields: A dictionary mapping obs_column to registry_field.
        feature_field: The registry field to validate variables index against.
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        The registered Artifact.
    """
    if isinstance(data, ad.AnnData):
        artifact = Artifact.from_anndata(data, description=description)
        artifact.n_observations = data.n_obs
    elif isinstance(data, pd.DataFrame):
        artifact = Artifact.from_df(data, description=description)
    else:
        raise ValueError("data must be a DataFrame or AnnData object")
    artifact.save()

    feature_kwargs: Dict = {}
    organism = check_registry_organism(
        feature_field.field.model, kwargs.pop("organism", None)
    )
    if organism is not None:
        feature_kwargs["organism"] = organism

    if isinstance(data, ad.AnnData):
        artifact.features.add_from_anndata(var_field=feature_field, **feature_kwargs)
    else:
        artifact.features.add_from_df(field=feature_field, **feature_kwargs)

    features = Feature.lookup().dict()
    for feature_name, field in fields.items():
        feature = features.get(feature_name)
        registry = field.field.model
        filter_kwargs = kwargs.copy()
        organism = check_registry_organism(registry, organism)
        if organism is not None:
            filter_kwargs["organism"] = organism
        df = data.obs if isinstance(data, ad.AnnData) else data
        labels = registry.from_values(df[feature_name], field=field, **filter_kwargs)
        artifact.labels.add(labels, feature)

    slug = ln_setup.settings.instance.slug
    logger.success(f"registered artifact in {colors.italic(slug)}")
    if ln_setup.settings.instance.is_remote:
        logger.info(f"🔗 https://lamin.ai/{slug}/artifact/{artifact.uid}")

    return artifact


def update_registry(
    values: List[str],
    field: FieldAttr,
    feature_name: str,
    using: Optional[str] = None,
    validated_only: bool = True,
    kwargs: Optional[Dict] = None,
    df: Optional[pd.DataFrame] = None,
) -> None:
    """Register features or labels records in the default instance from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        feature_name: The name of the feature to register.
        using: The name of the instance from which to transfer labels (if applicable).
        validated_only: If True, only register validated labels.
        kwargs: Additional keyword arguments to pass to the registry model.
        df: A DataFrame to register labels from.
    """
    from lamindb._save import save as ln_save
    from lamindb.core._settings import settings

    filter_kwargs = {} if kwargs is None else kwargs.copy()
    registry = field.field.model
    if registry == ULabel:
        validated_only = False

    organism = check_registry_organism(registry, filter_kwargs.pop("organism", None))
    if organism is not None:
        filter_kwargs["organism"] = organism

    verbosity = settings.verbosity
    try:
        settings.verbosity = "error"
        inspect_result_current = standardize_and_inspect(
            values=values, field=field, registry=registry, **filter_kwargs
        )
        if not inspect_result_current.non_validated:
            settings.verbosity = verbosity
            return

        labels_registered: Dict = {"from public": [], "without reference": []}

        (
            labels_registered[f"from {using}"],
            non_validated_labels,
        ) = update_registry_from_using_instance(
            inspect_result_current.non_validated,
            field=field,
            using=using,
            kwargs=filter_kwargs,
        )

        public_records = (
            registry.from_values(non_validated_labels, field=field, **filter_kwargs)
            if non_validated_labels
            else []
        )
        ln_save(public_records)
        labels_registered["from public"] = [
            getattr(r, field.field.name) for r in public_records
        ]
        labels_registered["without reference"] = [
            i for i in non_validated_labels if i not in labels_registered["from public"]
        ]

        if not validated_only:
            non_validated_records = []
            if df is not None and registry == Feature:
                non_validated_records = Feature.from_df(df)
            else:
                if "organism" in filter_kwargs:
                    filter_kwargs["organism"] = _register_organism(name=organism)
                for value in labels_registered["without reference"]:
                    filter_kwargs[field.field.name] = value
                    if registry == Feature:
                        filter_kwargs["type"] = "category"
                    non_validated_records.append(registry(**filter_kwargs))
            ln_save(non_validated_records)

        if registry == ULabel and field.field.name == "name":
            register_ulabels_with_parent(values, field=field, feature_name=feature_name)
    finally:
        settings.verbosity = verbosity

    log_registered_labels(
        labels_registered,
        feature_name=feature_name,
        model_field=f"{registry.__name__}.{field.field.name}",
        validated_only=validated_only,
    )


def log_registered_labels(
    labels_registered: Dict,
    feature_name: str,
    model_field: str,
    validated_only: bool = True,
) -> None:
    """Log the registered labels."""
    labels_type = "features" if feature_name == "feature" else "labels"
    model_field = colors.italic(model_field)
    for key, labels in labels_registered.items():
        if not labels:
            continue

        if key == "without reference" and validated_only:
            msg = colors.yellow(
                f"{len(labels)} non-validated {labels_type} are not registered with {model_field}: {labels}!"
            )
            lookup_print = f".lookup().['{feature_name}']"
            msg += f"\n      → to lookup categories, use {lookup_print}"
            msg += (
                f"\n      → to register, run {colors.yellow('register_features(validated_only=False)')}"
                if labels_type == "features"
                else f"\n      → to register, set {colors.yellow('validated_only=False')}"
            )
            logger.warning(msg)
        else:
            key = "" if key == "without reference" else f"{colors.green(key)} "
            logger.success(
                f"registered {len(labels)} {labels_type} {key}with {model_field}: {labels}"
            )


def register_ulabels_with_parent(
    values: List[str], field: FieldAttr, feature_name: str
) -> None:
    """Register a parent label for the given labels."""
    registry = field.field.model
    assert registry == ULabel
    all_records = registry.from_values(values, field=field)
    is_feature = registry.filter(name=f"is_{feature_name}").one_or_none()
    if is_feature is None:
        is_feature = registry(name=f"is_{feature_name}")
        is_feature.save()
    is_feature.children.add(*all_records)


def update_registry_from_using_instance(
    values: List[str],
    field: FieldAttr,
    using: Optional[str] = None,
    kwargs: Optional[Dict] = None,
) -> Tuple[List[str], List[str]]:
    """Register features or labels records from the using instance.

    Args:
        values: A list of values to be registered as labels.
        field: The FieldAttr object representing the field for which labels are being registered.
        using: The name of the instance from which to transfer labels (if applicable).
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        A tuple containing the list of registered labels and the list of non-registered labels.
    """
    kwargs = kwargs or {}
    labels_registered = []
    not_registered = values

    if using is not None and using != "default":
        registry = field.field.model
        registry_using = get_registry_instance(registry, using)
        inspect_result_using = standardize_and_inspect(
            values=values, field=field, registry=registry_using, **kwargs
        )
        labels_using = registry_using.filter(
            **{f"{field.field.name}__in": inspect_result_using.validated}
        ).all()
        for label_using in labels_using:
            label_using.save()
            labels_registered.append(getattr(label_using, field.field.name))
        not_registered = inspect_result_using.non_validated

    return labels_registered, not_registered


def _register_organism(name: str):
    """Register an organism record."""
    import bionty as bt

    organism = bt.Organism.filter(name=name).one_or_none()
    if organism is None:
        organism = bt.Organism.from_public(name=name)
        if organism is None:
            raise ValueError(
                f"Organism '{name}' not found\n"
                f"      → please register it: bt.Organism(name='{name}').save()"
            )
        organism.save()
    return organism
