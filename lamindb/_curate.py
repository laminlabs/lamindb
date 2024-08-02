from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lnschema_core import (
    Artifact,
    Feature,
    Record,
    Run,
    ULabel,
)

from .core.exceptions import ValidationError

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from lnschema_core.types import FieldAttr
    from mudata import MuData


class CurateLookup:
    """Lookup categories from the reference instance."""

    def __init__(
        self,
        categoricals: dict[str, FieldAttr],
        slots: dict[str, FieldAttr] = None,
        using: str | None = None,
    ) -> None:
        if slots is None:
            slots = {}
        self._fields = {**categoricals, **slots}
        self._using = None if using == "default" else using
        self._using_name = self._using or ln_setup.settings.instance.slug
        debug_message = f"Lookup objects from the " f"{colors.italic(self._using_name)}"
        logger.debug(debug_message)

    def __getattr__(self, name):
        if name in self._fields:
            registry = self._fields[name].field.model
            if self._using == "public":
                return registry.public().lookup()
            else:
                return get_registry_instance(registry, self._using).lookup()
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

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
            getattr_keys = "\n ".join(
                [f".{key}" for key in self._fields if key.isidentifier()]
            )
            getitem_keys = "\n ".join(
                [str([key]) for key in self._fields if not key.isidentifier()]
            )
            return (
                f"Lookup objects from the {colors.italic(self._using_name)}:\n "
                f"{colors.green(getattr_keys)}\n "
                f"{colors.green(getitem_keys)}\n\n"
                "Example:\n    → categories = validator.lookup().cell_type\n"
                "    → categories.alveolar_type_1_fibroblast_cell"
            )
        else:  # pragma: no cover
            return colors.warning("No fields are found!")


class DataFrameCurator:
    """Annotation flow for a DataFrame object.

    Args:
        df: The DataFrame object to curate.
        columns: The field attribute for the feature column.
        categoricals: A dictionary mapping column names to registry_field.
        using: The reference instance containing registries to validate against.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping column names to Source records.

    Examples:
        >>> import bionty as bt
        >>> curate = ln.Curate.from_df(
                df,
                categoricals={"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name}
            )
    """

    def __init__(
        self,
        df: pd.DataFrame,
        columns: FieldAttr = Feature.name,
        categoricals: dict[str, FieldAttr] | None = None,
        using: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> None:
        from lamindb.core._settings import settings

        self._df = df
        self._fields = categoricals or {}
        self._columns_field = columns
        self._using = using
        settings.verbosity = verbosity
        self._artifact = None
        self._collection = None
        self._validated = False
        self._kwargs = {"organism": organism} if organism else {}
        if sources is None:
            sources = {}
        self._sources = sources
        self._save_columns()

    @property
    def fields(self) -> dict:
        """Return the columns fields to validate against."""
        return self._fields

    def lookup(self, using: str | None = None) -> CurateLookup:
        """Lookup categories.

        Args:
            using: The instance where the lookup is performed.
                if None (default), the lookup is performed on the instance specified in "using" parameter of the validator.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._fields,
            slots={"columns": self._columns_field},
            using=using or self._using,
        )

    def _save_columns(self, validated_only: bool = True, **kwargs) -> None:
        """Save column name records."""
        missing_columns = set(self.fields.keys()) - set(self._df.columns)
        if missing_columns:
            raise ValueError(
                f"Columns {missing_columns} are not found in the data object!"
            )

        # Always save features specified as the fields keys
        update_registry(
            values=list(self.fields.keys()),
            field=self._columns_field,
            key="columns",
            save_function="add_new_from_columns",
            using=self._using,
            validated_only=False,
            source=self._sources.get("columns"),
            **kwargs,
        )

        # Save the rest of the columns based on validated_only
        additional_columns = set(self._df.columns) - set(self.fields.keys())
        if additional_columns:
            update_registry(
                values=list(additional_columns),
                field=self._columns_field,
                key="columns",
                save_function="add_new_from_columns",
                using=self._using,
                validated_only=validated_only,
                df=self._df,  # Get the Feature type from df
                source=self._sources.get("columns"),
                **kwargs,
            )

    def add_validated_from(self, key: str, organism: str | None = None):
        """Add validated categories.

        Args:
            key: The key referencing the slot in the DataFrame.
            organism: The organism name.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._update_registry(key, validated_only=True, **self._kwargs)

    def add_new_from(self, key: str, organism: str | None = None, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        self._kwargs.update({"organism": organism} if organism else {})
        self._update_registry(key, validated_only=False, **self._kwargs, **kwargs)

    def add_new_from_columns(self, organism: str | None = None, **kwargs):
        """Add validated & new column names to its registry.

        Args:
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._save_columns(validated_only=False, **self._kwargs, **kwargs)

    def _update_registry(self, categorical: str, validated_only: bool = True, **kwargs):
        if categorical == "all":
            self._update_registry_all(validated_only=validated_only, **kwargs)
        elif categorical == "columns":
            self._save_columns(validated_only=validated_only, **kwargs)
        else:
            if categorical not in self.fields:
                raise ValueError(f"Feature {categorical} is not part of the fields!")
            update_registry(
                values=self._df[categorical].unique().tolist(),
                field=self.fields[categorical],
                key=categorical,
                using=self._using,
                validated_only=validated_only,
                sources=self._sources.get(categorical),
                **kwargs,
            )

    def _update_registry_all(self, validated_only: bool = True, **kwargs):
        """Save labels for all features."""
        for name in self.fields.keys():
            logger.info(f"saving labels for '{name}'")
            self._update_registry(name, validated_only=validated_only, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate variables and categorical observations.

        Returns:
            Whether the DataFrame is validated.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._validated = validate_categories_in_df(
            self._df,
            fields=self.fields,
            using=self._using,
            sources=self._sources,
            **self._kwargs,
        )
        return self._validated

    def save_artifact(self, description: str | None = None, **kwargs) -> Artifact:
        """Save the validated DataFrame and metadata.

        Args:
            description: Description of the DataFrame object.
            **kwargs: Object level metadata.

        Returns:
            A saved artifact record.
        """
        from lamindb.core._settings import settings

        if not self._validated:
            raise ValidationError(
                f"Data object is not validated, please run {colors.yellow('validate()')}!"
            )

        # Make sure all labels are saved in the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            # save all validated records to the current instance
            self.add_validated_from("all")

            self._artifact = save_artifact(
                self._df,
                description=description,
                fields=self.fields,
                columns_field=self._columns_field,
                **kwargs,
                **self._kwargs,
            )
        finally:
            settings.verbosity = verbosity

        return self._artifact

    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't save any outputs."""
        from lamindb.core._run_context import run_context

        if run_context.transform is not None:
            Run.filter(transform=run_context.transform, output_artifacts=None).exclude(
                uid=run_context.run.uid
            ).delete()


class AnnDataCurator(DataFrameCurator):
    """Annotation flow for ``AnnData``.

    Args:
        data: The AnnData object or an AnnData-like path.
        var_index: The registry field for mapping the ``.var`` index.
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
        using: A reference LaminDB instance.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping ``.obs.columns`` to Source records.

    Examples:
        >>> import bionty as bt
        >>> curate = ln.Curate.from_anndata(
                adata,
                var_index=bt.Gene.ensembl_gene_id,
                categoricals={"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name},
                organism="human",
            )
    """

    def __init__(
        self,
        data: ad.AnnData | UPathStr,
        var_index: FieldAttr,
        categoricals: dict[str, FieldAttr] | None = None,
        using: str = "default",
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> None:
        from lamindb_setup.core import upath

        from ._artifact import data_is_anndata

        if sources is None:
            sources = {}
        if not data_is_anndata(data):
            raise ValueError(
                "data has to be an AnnData object or a path to AnnData-like"
            )
        if isinstance(data, ad.AnnData):
            self._adata = data
        else:  # pragma: no cover
            from lamindb.core.storage._backed_access import backed_access

            self._adata = backed_access(upath.create_path(data))

        self._data = data
        self._var_field = var_index
        super().__init__(
            df=self._adata.obs,
            categoricals=categoricals,
            using=using,
            verbosity=verbosity,
            organism=organism,
            sources=sources,
        )
        self._obs_fields = categoricals
        self._save_from_var_index(validated_only=True, **self._kwargs)

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(self, using: str | None = None) -> CurateLookup:
        """Lookup categories.

        Args:
            using: The instance where the lookup is performed.
                if None (default), the lookup is performed on the instance specified in "using" parameter of the validator.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, "var_index": self._var_field},
            using=using or self._using,
        )

    def _save_from_var_index(
        self, validated_only: bool = True, organism: str | None = None
    ):
        """Save variable records."""
        update_registry(
            values=self._adata.var.index,
            field=self.var_index,
            key="var_index",
            save_function="add_new_from_var_index",
            using=self._using,
            validated_only=validated_only,
            organism=organism,
            source=self._sources.get("var_index"),
        )

    def add_new_from_var_index(self, organism: str | None = None, **kwargs):
        """Update variable records.

        Args:
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._save_from_var_index(validated_only=False, **self._kwargs, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate categories.

        Args:
            organism: The organism name.

        Returns:
            Whether the AnnData object is validated.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        if self._using is not None and self._using != "default":
            logger.important(
                f"validating metadata using registries of instance {colors.italic(self._using)}"
            )
        validated_var = validate_categories(
            self._adata.var.index,
            field=self._var_field,
            key="var_index",
            using=self._using,
            **self._kwargs,
        )
        validated_obs = validate_categories_in_df(
            self._adata.obs,
            fields=self.categoricals,
            using=self._using,
            sources=self._sources,
            **self._kwargs,
        )
        self._validated = validated_var and validated_obs
        return self._validated

    def save_artifact(self, description: str | None = None, **kwargs) -> Artifact:
        """Save the validated ``AnnData`` and metadata.

        Args:
            description: Description of the ``AnnData`` object.
            **kwargs: Object level metadata.

        Returns:
            A saved artifact record.
        """
        if not self._validated:
            raise ValidationError(
                f"Data object is not validated, please run {colors.yellow('validate()')}!"
            )

        self._artifact = save_artifact(
            self._data,
            adata=self._adata,
            description=description,
            columns_field=self.var_index,
            fields=self.categoricals,
            **self._kwargs,
            **kwargs,
        )
        return self._artifact


class MuDataCurator:
    """Annotation flow for a ``MuData`` object.

    Args:
        mdata: The MuData object to curate.
        var_index: The registry field for mapping the ``.var`` index for each modality.
            For example:
            ``{"modality_1": bt.Gene.ensembl_gene_id, "modality_2": ln.CellMarker.name}``
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
            Use modality keys to specify categoricals for MuData slots such as `"rna:cell_type": bt.CellType.name"`.
        using: A reference LaminDB instance.
        verbosity: The verbosity level.
        organism: The organism name.

    Examples:
        >>> import bionty as bt
        >>> curate = ln.Curate.from_mudata(
                mdata,
                var_index={"rna": bt.Gene.ensembl_gene_id, "adt": ln.CellMarker.name},
                categoricals={"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ln.ULabel.name},
                organism="human",
            )
    """

    def __init__(
        self,
        mdata: MuData,
        var_index: dict[str, dict[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        using: str = "default",
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> None:
        if sources is None:
            sources = {}
        self._sources = sources
        self._mdata = mdata
        self._kwargs = {"organism": organism} if organism else {}
        self._var_fields = var_index
        self._verify_modality(self._var_fields.keys())
        self._obs_fields = self._parse_categoricals(categoricals)
        self._modalities = set(self._var_fields.keys()) | set(self._obs_fields.keys())
        self._using = using
        self._verbosity = verbosity
        self._df_annotators = {
            modality: DataFrameCurator(
                df=mdata[modality].obs if modality != "obs" else mdata.obs,
                categoricals=self._obs_fields.get(modality, {}),
                using=using,
                verbosity=verbosity,
                sources=self._sources.get(modality),
                **self._kwargs,
            )
            for modality in self._modalities
        }
        for modality in self._var_fields.keys():
            self._save_from_var_index_modality(
                modality=modality, validated_only=True, **self._kwargs
            )

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_fields

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def _verify_modality(self, modalities: Iterable[str]):
        """Verify the modality exists."""
        for modality in modalities:
            if modality not in self._mdata.mod.keys():
                raise ValueError(f"modality '{modality}' does not exist!")

    def _save_from_var_index_modality(
        self, modality: str, validated_only: bool = True, **kwargs
    ):
        """Save variable records."""
        update_registry(
            values=self._mdata[modality].var.index,
            field=self._var_fields[modality],
            key="var_index",
            save_function="add_new_from_var_index",
            using=self._using,
            validated_only=validated_only,
            dtype="number",
            **kwargs,
        )

    def _parse_categoricals(self, categoricals: dict[str, FieldAttr]) -> dict:
        """Parse the categorical fields."""
        prefixes = {f"{k}:" for k in self._mdata.mod.keys()}
        obs_fields: dict[str, dict[str, FieldAttr]] = {}
        for k, v in categoricals.items():
            if k not in self._mdata.obs.columns:
                raise ValueError(f"column '{k}' does not exist in mdata.obs!")
            if any(k.startswith(prefix) for prefix in prefixes):
                modality, col = k.split(":")[0], k.split(":")[1]
                if modality not in obs_fields.keys():
                    obs_fields[modality] = {}
                obs_fields[modality][col] = v
            else:
                if "obs" not in obs_fields.keys():
                    obs_fields["obs"] = {}
                obs_fields["obs"][k] = v
        return obs_fields

    def lookup(self, using: str | None = None) -> CurateLookup:
        """Lookup categories.

        Args:
            using: The instance where the lookup is performed.
                if None (default), the lookup is performed on the instance specified in "using" parameter of the validator.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={
                **self._obs_fields,
                **{f"{k}_var_index": v for k, v in self._var_fields.items()},
            },
            using=using or self._using,
        )

    def add_new_from_columns(
        self,
        modality: str,
        column_names: list[str] | None = None,
        organism: str | None = None,
        **kwargs,
    ):
        """Update columns records.

        Args:
            modality: The modality name.
            column_names: The column names to save.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        update_registry(
            values=column_names or self._mdata[modality].obs.columns,
            field=Feature.name,
            key=f"{modality} obs columns",
            using=self._using,
            validated_only=False,
            df=self._mdata[modality].obs,
            **self._kwargs,
            **kwargs,
        )

    def add_new_from_var_index(
        self, modality: str, organism: str | None = None, **kwargs
    ):
        """Update variable records.

        Args:
            modality: The modality name.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._save_from_var_index_modality(
            modality=modality, validated_only=False, **self._kwargs, **kwargs
        )

    def add_validated_from(
        self, key: str, modality: str | None = None, organism: str | None = None
    ):
        """Add validated categories.

        Args:
            key: The key referencing the slot in the DataFrame.
            modality: The modality name.
            organism: The organism name.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        modality = modality or "obs"
        if modality in self._df_annotators:
            df_annotator = self._df_annotators[modality]
            df_annotator.add_validated_from(key=key, **self._kwargs)

    def add_new_from(
        self,
        key: str,
        modality: str | None = None,
        organism: str | None = None,
        **kwargs,
    ):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame.
            modality: The modality name.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to the registry model.
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        self._kwargs.update({"organism": organism} if organism else {})
        modality = modality or "obs"
        if modality in self._df_annotators:
            df_annotator = self._df_annotators[modality]
            df_annotator.add_new_from(key=key, **self._kwargs, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate categories."""
        self._kwargs.update({"organism": organism} if organism else {})
        if self._using is not None and self._using != "default":
            logger.important(
                f"validating metadata using registries of instance {colors.italic(self._using)}"
            )
        validated_var = True
        for modality, var_field in self._var_fields.items():
            validated_var &= validate_categories(
                self._mdata[modality].var.index,
                field=var_field,
                key=f"{modality}_var_index",
                using=self._using,
                **self._kwargs,
            )
        validated_obs = True
        for modality, fields in self._obs_fields.items():
            if modality == "obs":
                obs = self._mdata.obs
            else:
                obs = self._mdata[modality].obs
            validated_obs &= validate_categories_in_df(
                obs,
                fields=fields,
                using=self._using,
                sources=self._sources.get(modality),
                **self._kwargs,
            )
        self._validated = validated_var and validated_obs
        return self._validated

    def save_artifact(self, description: str | None = None, **kwargs) -> Artifact:
        """Save the validated ``MuData`` and metadata.

        Args:
            description: Description of the ``MuData`` object.
            **kwargs: Object level metadata.

        Returns:
            A saved artifact record.
        """
        if not self._validated:
            raise ValidationError("Please run `validate()` first!")

        self._artifact = save_artifact(
            self._mdata,
            description=description,
            columns_field=self.var_index,
            fields=self.categoricals,
            **self._kwargs,
            **kwargs,
        )
        return self._artifact


class Curate:
    """Annotation flow."""

    @classmethod
    @doc_args(DataFrameCurator.__doc__)
    def from_df(
        cls,
        df: pd.DataFrame,
        categoricals: dict[str, FieldAttr] | None = None,
        columns: FieldAttr = Feature.name,
        using: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
    ) -> DataFrameCurator:
        """{}"""  # noqa: D415
        return DataFrameCurator(
            df=df,
            categoricals=categoricals,
            columns=columns,
            using=using,
            verbosity=verbosity,
            organism=organism,
        )

    @classmethod
    @doc_args(AnnDataCurator.__doc__)
    def from_anndata(
        cls,
        data: ad.AnnData | UPathStr,
        var_index: FieldAttr,
        categoricals: dict[str, FieldAttr] | None = None,
        using: str = "default",
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> AnnDataCurator:
        """{}"""  # noqa: D415
        return AnnDataCurator(
            data=data,
            var_index=var_index,
            categoricals=categoricals,
            using=using,
            verbosity=verbosity,
            organism=organism,
            sources=sources,
        )

    @classmethod
    @doc_args(MuDataCurator.__doc__)
    def from_mudata(
        cls,
        mdata: MuData,
        var_index: dict[str, dict[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        using: str = "default",
        verbosity: str = "hint",
        organism: str | None = None,
    ) -> MuDataCurator:
        """{}"""  # noqa: D415
        return MuDataCurator(
            mdata=mdata,
            var_index=var_index,
            categoricals=categoricals,
            using=using,
            verbosity=verbosity,
            organism=organism,
        )


def get_registry_instance(registry: Record, using: str | None = None) -> Record:
    """Get a registry instance using a specific instance."""
    if using is not None and using != "default":
        return registry.using(using)
    return registry


def standardize_and_inspect(
    values: Iterable[str], field: FieldAttr, registry: Record, **kwargs
):
    """Standardize and inspect values using a registry."""
    if hasattr(registry, "standardize") and hasattr(
        registry,
        "synonyms",  # https://github.com/laminlabs/lamindb/issues/1685
    ):
        values = registry.standardize(values, field=field, mute=True, **kwargs)
    return registry.inspect(values, field=field, mute=True, **kwargs)


def check_registry_organism(registry: Record, organism: str | None = None) -> dict:
    """Check if a registry needs an organism and return the organism name."""
    if hasattr(registry, "organism_id"):
        import bionty as bt

        if organism is None and bt.settings.organism is None:
            raise ValueError(
                f"{registry.__name__} registry requires an organism!\n"
                "      → please pass an organism name via organism="
            )
        return {"organism": organism or bt.settings.organism.name}
    return {}


def validate_categories(
    values: Iterable[str],
    field: FieldAttr,
    key: str,
    using: str | None = None,
    organism: str | None = None,
    source: Record | None = None,
) -> bool:
    """Validate ontology terms in a pandas series using LaminDB registries."""
    from lamindb._from_values import _print_values
    from lamindb.core._settings import settings

    model_field = f"{field.field.model.__name__}.{field.field.name}"

    def _log_mapping_info():
        logger.indent = ""
        logger.info(f"mapping {colors.italic(key)} on {colors.italic(model_field)}")
        logger.indent = "   "

    registry = field.field.model
    filter_kwargs = check_registry_organism(registry, organism)
    filter_kwargs.update({"source": source} if source else {})

    # Inspect the default instance
    inspect_result = standardize_and_inspect(
        values=values, field=field, registry=registry, **filter_kwargs
    )
    non_validated = inspect_result.non_validated

    values_validated = []
    if using is not None and using != "default" and non_validated:
        registry = get_registry_instance(registry, using)
        # Inspect the using instance
        inspect_result = standardize_and_inspect(
            values=non_validated, field=field, registry=registry, **filter_kwargs
        )
        non_validated = inspect_result.non_validated
        values_validated += inspect_result.validated

    # Inspect from public (bionty only)
    if hasattr(registry, "public"):
        verbosity = settings.verbosity
        try:
            settings.verbosity = "error"
            public_records = registry.from_values(
                non_validated, field=field, **filter_kwargs
            )
            values_validated += [getattr(r, field.field.name) for r in public_records]
        finally:
            settings.verbosity = verbosity

    validated_hint_print = f".add_validated_from('{key}')"
    n_validated = len(values_validated)
    if n_validated > 0:
        _log_mapping_info()
        logger.warning(
            f"found {colors.yellow(f'{n_validated} terms')} validated terms: "
            f"{colors.yellow(values_validated)}\n      → save terms via "
            f"{colors.yellow(validated_hint_print)}"
        )

    non_validated_hint_print = f".add_new_from('{key}')"
    non_validated = [i for i in non_validated if i not in values_validated]
    n_non_validated = len(non_validated)
    if n_non_validated == 0:
        logger.indent = ""
        logger.success(f"{key} is validated against {colors.italic(model_field)}")
        return True
    else:
        are = "are" if n_non_validated > 1 else "is"
        print_values = _print_values(non_validated)
        warning_message = (
            f"{colors.yellow(f'{n_non_validated} terms')} {are} not validated: "
            f"{colors.yellow(print_values)}\n      → save terms via "
            f"{colors.yellow(non_validated_hint_print)}"
        )
        if logger.indent == "":
            _log_mapping_info()
        logger.warning(warning_message)
        logger.indent = ""
        return False


def validate_categories_in_df(
    df: pd.DataFrame,
    fields: dict[str, FieldAttr],
    using: str | None = None,
    sources: dict[str, Record] = None,
    **kwargs,
) -> bool:
    """Validate categories in DataFrame columns using LaminDB registries."""
    if sources is None:
        sources = {}
    validated = True
    for key, field in fields.items():
        validated &= validate_categories(
            df[key],
            field=field,
            key=key,
            using=using,
            source=sources.get(key),
            **kwargs,
        )
    return validated


def save_artifact(
    data: pd.DataFrame | ad.AnnData | MuData,
    fields: dict[str, FieldAttr] | dict[str, dict[str, FieldAttr]],
    columns_field: FieldAttr | dict[str, FieldAttr],
    description: str | None = None,
    organism: str | None = None,
    adata: ad.AnnData | None = None,
    **kwargs,
) -> Artifact:
    """Save all metadata with an Artifact.

    Args:
        data: The DataFrame or AnnData object to save.
        description: A description of the artifact.
        fields: A dictionary mapping obs_column to registry_field.
        columns_field: The registry field to validate variables index against.
        organism: The organism name.
        adata: The AnnData object to save, must be provided if data is a path.
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        The saved Artifact.
    """
    from ._artifact import data_is_anndata

    artifact = None
    if data_is_anndata(data):
        assert adata is not None  # noqa: S101
        artifact = Artifact.from_anndata(data, description=description, **kwargs)
        artifact.n_observations = adata.shape[0]
        data = adata

    elif isinstance(data, pd.DataFrame):
        artifact = Artifact.from_df(data, description=description, **kwargs)
    else:
        try:
            from mudata import MuData

            if isinstance(data, MuData):
                artifact = Artifact.from_mudata(data, description=description, **kwargs)
                artifact.n_observations = data.n_obs
        except ImportError:
            pass
    if artifact is None:
        raise ValueError("data must be a DataFrame, AnnData or MuData object.")
    artifact.save()

    feature_kwargs = check_registry_organism(
        (
            list(columns_field.values())[0].field.model
            if isinstance(columns_field, dict)
            else columns_field.field.model
        ),
        organism,
    )

    if artifact._accessor == "DataFrame":
        artifact.features._add_set_from_df(field=columns_field, **feature_kwargs)
    elif artifact._accessor == "AnnData":
        artifact.features._add_set_from_anndata(
            var_field=columns_field, **feature_kwargs
        )
    elif artifact._accessor == "MuData":
        artifact.features._add_set_from_mudata(
            var_fields=columns_field, **feature_kwargs
        )
    else:
        raise NotImplementedError

    def _add_labels(data, artifact: Artifact, fields: dict[str, FieldAttr]):
        features = Feature.lookup().dict()
        for key, field in fields.items():
            feature = features.get(key)
            registry = field.field.model
            filter_kwargs = check_registry_organism(registry, organism)
            df = data if isinstance(data, pd.DataFrame) else data.obs
            labels = registry.from_values(df[key], field=field, **filter_kwargs)
            artifact.labels.add(labels, feature)

    if artifact._accessor == "MuData":
        for modality, modality_fields in fields.items():
            if modality == "obs":
                _add_labels(data, artifact, modality_fields)
            else:
                _add_labels(data[modality], artifact, modality_fields)
    else:
        _add_labels(data, artifact, fields)

    slug = ln_setup.settings.instance.slug
    if ln_setup.settings.instance.is_remote:  # pragma: no cover
        logger.important(f"go to https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact


def update_registry(
    values: list[str],
    field: FieldAttr,
    key: str,
    save_function: str = "add_new_from",
    using: str | None = None,
    validated_only: bool = True,
    df: pd.DataFrame | None = None,
    organism: str | None = None,
    dtype: str | None = None,
    source: Record | None = None,
    **kwargs,
) -> list[Record]:
    """Save features or labels records in the default instance from the using instance.

    Args:
        values: A list of values to be saved as labels.
        field: The FieldAttr object representing the field for which labels are being saved.
        key: The name of the feature to save.
        save_function: The name of the function to save the labels.
        using: The name of the instance from which to transfer labels (if applicable).
        validated_only: If True, only save validated labels.
        df: A DataFrame to save labels from.
        organism: The organism name.
        dtype: The type of the feature.
        source: The source record.
        kwargs: Additional keyword arguments to pass to the registry model to create new records.
    """
    from lamindb._save import save as ln_save
    from lamindb.core._settings import settings

    registry = field.field.model
    filter_kwargs = check_registry_organism(registry, organism)
    filter_kwargs.update({"source": source} if source else {})

    verbosity = settings.verbosity
    try:
        settings.verbosity = "error"
        inspect_result_current = standardize_and_inspect(
            values=values, field=field, registry=registry, **filter_kwargs
        )
        if not inspect_result_current.non_validated:
            all_labels = registry.from_values(
                inspect_result_current.validated, field=field, **filter_kwargs
            )
            settings.verbosity = verbosity
            return all_labels

        labels_saved: dict = {"from public": [], "without reference": []}

        (
            labels_saved[f"from {using}"],
            non_validated_labels,
        ) = update_registry_from_using_instance(
            inspect_result_current.non_validated,
            field=field,
            using=using,
            **filter_kwargs,
        )

        public_records = (
            registry.from_values(non_validated_labels, field=field, **filter_kwargs)
            if non_validated_labels
            else []
        )
        # here we check to only save the public records if they are from the specified source
        # TODO: this if shouldn't be needed
        if source:
            public_records = [r for r in public_records if r.source == source]
        ln_save(public_records)
        labels_saved["from public"] = [
            getattr(r, field.field.name) for r in public_records
        ]
        labels_saved["without reference"] = [
            i for i in non_validated_labels if i not in labels_saved["from public"]
        ]

        if not validated_only:
            non_validated_records = []
            if df is not None and registry == Feature:
                non_validated_records = Feature.from_df(df)
            else:
                if "organism" in filter_kwargs:
                    filter_kwargs["organism"] = _save_organism(name=organism)
                init_kwargs = {}
                for value in labels_saved["without reference"]:
                    init_kwargs[field.field.name] = value
                    if registry == Feature:
                        init_kwargs["dtype"] = "cat" if dtype is None else dtype
                    non_validated_records.append(
                        registry(
                            **init_kwargs,
                            **{k: v for k, v in filter_kwargs.items() if k != "source"},
                            **{k: v for k, v in kwargs.items() if k != "sources"},
                        )
                    )
            ln_save(non_validated_records)

        if registry == ULabel and field.field.name == "name":
            save_ulabels_with_parent(values, field=field, key=key)

        # get all records
        all_labels = registry.from_values(
            inspect_result_current.validated + inspect_result_current.non_validated,
            field=field,
            **filter_kwargs,
        )
    finally:
        settings.verbosity = verbosity

    log_saved_labels(
        labels_saved,
        key=key,
        save_function=save_function,
        model_field=f"{registry.__name__}.{field.field.name}",
        validated_only=validated_only,
    )

    return all_labels


def log_saved_labels(
    labels_saved: dict,
    key: str,
    save_function: str,
    model_field: str,
    validated_only: bool = True,
) -> None:
    """Log the saved labels."""
    from ._from_values import _print_values

    model_field = colors.italic(model_field)
    for k, labels in labels_saved.items():
        if not labels:
            continue

        if k == "without reference" and validated_only:
            msg = colors.yellow(
                f"{len(labels)} non-validated categories are not saved in {model_field}: {labels}!"
            )
            lookup_print = (
                f"lookup().{key}" if key.isidentifier() else f".lookup()['{key}']"
            )

            hint = f".add_new_from('{key}')"
            msg += f"\n      → to lookup categories, use {lookup_print}"
            msg += (
                f"\n      → to save, run {colors.yellow(hint)}"
                if save_function == "add_new_from"
                else f"\n      → to save, run {colors.yellow(save_function)}"
            )
            logger.warning(msg)
        else:
            k = "" if k == "without reference" else f"{colors.green(k)} "
            # the term "transferred" stresses that this is always in the context of transferring
            # labels from a public ontology or a different instance to the present instance
            s = "s" if len(labels) > 1 else ""
            logger.success(
                f"added {len(labels)} record{s} {k}with {model_field} for {colors.italic(key)}: {_print_values(labels)}"
            )


def save_ulabels_with_parent(values: list[str], field: FieldAttr, key: str) -> None:
    """Save a parent label for the given labels."""
    registry = field.field.model
    assert registry == ULabel  # noqa: S101
    all_records = registry.from_values(values, field=field)
    is_feature = registry.filter(name=f"is_{key}").one_or_none()
    if is_feature is None:
        is_feature = registry(name=f"is_{key}")
        is_feature.save()
    is_feature.children.add(*all_records)


def update_registry_from_using_instance(
    values: list[str],
    field: FieldAttr,
    using: str | None = None,
    **kwargs,
) -> tuple[list[str], list[str]]:
    """Save features or labels records from the using instance.

    Args:
        values: A list of values to be saved as labels.
        field: The FieldAttr object representing the field for which labels are being saved.
        using: The name of the instance from which to transfer labels (if applicable).
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        A tuple containing the list of saved labels and the list of non-saved labels.
    """
    labels_saved = []
    not_saved = values

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
            labels_saved.append(getattr(label_using, field.field.name))
        not_saved = inspect_result_using.non_validated

    return labels_saved, not_saved


def _save_organism(name: str):  # pragma: no cover
    """Save an organism record."""
    import bionty as bt

    organism = bt.Organism.filter(name=name).one_or_none()
    if organism is None:
        organism = bt.Organism.from_source(name=name)
        if organism is None:
            raise ValueError(
                f"Organism '{name}' not found\n"
                f"      → please save it: bt.Organism(name='{name}').save()"
            )
        organism.save()
    return organism
