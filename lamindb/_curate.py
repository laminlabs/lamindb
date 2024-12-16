from __future__ import annotations

import copy
import warnings
from itertools import chain
from typing import TYPE_CHECKING

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
import pyarrow as pa
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.upath import UPath
from lnschema_core import (
    Artifact,
    Feature,
    FeatureSet,
    Record,
    Run,
    ULabel,
)

from ._from_values import _format_values
from .core.exceptions import ValidationError

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any

    from lamindb_setup.core.types import UPathStr
    from lnschema_core.types import FieldAttr
    from mudata import MuData


class CurateLookup:
    """Lookup categories from the reference instance.

    Args:
        categoricals: A dictionary of categorical fields to lookup.
        slots: A dictionary of slot fields to lookup.
        using_key: The key of the instance to lookup from. Defaults to the
            current instance if not specified.
        public: Whether to lookup from the public instance. Defaults to False.

    Example:
        >>> curator = ln.Curator.from_df(...)
        >>> curator.lookup()["cell_type"].alveolar_type_1_fibroblast_cell
        <Category: alveolar_type_1_fibroblast_cell>

    """

    def __init__(
        self,
        categoricals: dict[str, FieldAttr],
        slots: dict[str, FieldAttr] = None,
        using_key: str | None = None,
        public: bool = False,
    ) -> None:
        slots = slots or {}
        self._fields = {**categoricals, **slots}
        self._using_key = None if using_key == "default" else using_key
        self._using_key_name = self._using_key or ln_setup.settings.instance.slug
        self._public = public
        debug_message = f"Lookup objects from {colors.italic(self._using_key_name)}"
        logger.debug(debug_message)

    def __getattr__(self, name):
        if name in self._fields:
            registry = self._fields[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public().lookup()
            else:
                return get_registry_instance(registry, self._using_key).lookup()
        raise AttributeError(
            f'"{self.__class__.__name__}" object has no attribute "{name}"'
        )

    def __getitem__(self, name):
        if name in self._fields:
            registry = self._fields[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public().lookup()
            else:
                return get_registry_instance(registry, self._using_key).lookup()
        raise AttributeError(
            f'"{self.__class__.__name__}" object has no attribute "{name}"'
        )

    def __repr__(self) -> str:
        if len(self._fields) > 0:
            getattr_keys = "\n ".join(
                [f".{key}" for key in self._fields if key.isidentifier()]
            )
            getitem_keys = "\n ".join(
                [str([key]) for key in self._fields if not key.isidentifier()]
            )
            ref = "public" if self._public else self._using_key_name
            return (
                f"Lookup objects from the {colors.italic(ref)}:\n "
                f"{colors.green(getattr_keys)}\n "
                f"{colors.green(getitem_keys)}\n"
                'Example:\n    → categories = curator.lookup()["cell_type"]\n'
                "    → categories.alveolar_type_1_fibroblast_cell\n\n"
                "To look up public ontologies, use .lookup(public=True)"
            )
        else:  # pragma: no cover
            return colors.warning("No fields are found!")


class BaseCurator:
    """Curate a dataset."""

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        import sys

        # Deprecated methods
        if "sphinx" not in sys.modules:
            if hasattr(cls, "_add_new_from_columns"):
                cls.add_new_from_columns = cls._add_new_from_columns

    def validate(self) -> bool:
        """Validate dataset.

        This method also registers the validated records in the current instance.

        Returns:
            Boolean indicating whether the dataset is validated.
        """
        pass  # pragma: no cover

    def standardize(self, key: str) -> None:
        """Replace synonyms with standardized values.

        Inplace modification of the dataset.

        Args:
            key: The name of the column to standardize.

        Returns:
            None
        """
        pass  # pragma: no cover

    def save_artifact(
        self,
        description: str | None = None,
        key: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Save the dataset as artifact.

        Args:
            description: A description of the DataFrame object.
            key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        pass  # pragma: no cover


class DataFrameCurator(BaseCurator):
    """Curation flow for a DataFrame object.

    See also :class:`~lamindb.Curator`.

    Args:
        df: The DataFrame object to curate.
        columns: The field attribute for the feature column.
        categoricals: A dictionary mapping column names to registry_field.
        using_key: The reference instance containing registries to validate against.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping column names to Source records.
        exclude: A dictionary mapping column names to values to exclude from validation.
            When specific :class:`~bionty.Source` instances are pinned and may lack default values (e.g., "unknown" or "na"),
            using the exclude parameter ensures they are not validated.

    Returns:
        A curator object.

    Examples:
        >>> import bionty as bt
        >>> curator = ln.Curator.from_df(
        ...     df,
        ...     categoricals={
        ...         "cell_type_ontology_id": bt.CellType.ontology_id,
        ...         "donor_id": ln.ULabel.name
        ...     }
        ... )
    """

    def __init__(
        self,
        df: pd.DataFrame,
        columns: FieldAttr = Feature.name,
        categoricals: dict[str, FieldAttr] | None = None,
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict | None = None,
        check_valid_keys: bool = True,
    ) -> None:
        from lamindb.core._settings import settings

        self._df = df
        self._fields = categoricals or {}
        self._columns_field = columns
        self._using_key = using_key
        # TODO: change verbosity back
        settings.verbosity = verbosity
        self._artifact = None
        self._collection = None
        self._validated = False
        self._kwargs = {"organism": organism} if organism else {}
        self._sources = sources or {}
        self._exclude = exclude or {}
        self._non_validated = None
        if check_valid_keys:
            self._check_valid_keys()
        self._save_columns()

    @property
    def non_validated(self) -> dict[str, list[str]]:
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return self._non_validated

    @property
    def fields(self) -> dict:
        """Return the columns fields to validate against."""
        return self._fields

    def lookup(
        self, using_key: str | None = None, public: bool = False
    ) -> CurateLookup:
        """Lookup categories.

        Args:
            using_key: The instance where the lookup is performed.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._fields,
            slots={"columns": self._columns_field},
            using_key=using_key or self._using_key,
            public=public,
        )

    def _check_valid_keys(self, extra: set | None = None) -> None:
        extra = extra or set()
        for name, d in {
            "categoricals": self._fields,
            "sources": self._sources,
            "exclude": self._exclude,
        }.items():
            if not isinstance(d, dict):
                raise TypeError(f"{name} must be a dictionary!")
            valid_keys = set(self._df.columns) | {"columns"} | extra
            nonval_keys = [key for key in d.keys() if key not in valid_keys]
            n = len(nonval_keys)
            s = "s" if n > 1 else ""
            are = "are" if n > 1 else "is"
            if len(nonval_keys) > 0:
                raise ValidationError(
                    f"key{s} passed to {name} {are} not present in columns: {colors.yellow(_format_values(nonval_keys))}"
                )

    def _save_columns(self, validated_only: bool = True) -> None:
        """Save column name records."""
        # Always save features specified as the fields keys
        update_registry(
            values=list(self.fields.keys()),
            field=self._columns_field,
            key="columns",
            using_key=self._using_key,
            validated_only=False,
            source=self._sources.get("columns"),
            exclude=self._exclude.get("columns"),
            **self._kwargs,  # type: ignore
        )

        # Save the rest of the columns based on validated_only
        additional_columns = set(self._df.columns) - set(self.fields.keys())
        if additional_columns:
            update_registry(
                values=list(additional_columns),
                field=self._columns_field,
                key="columns",
                using_key=self._using_key,
                validated_only=validated_only,
                df=self._df,  # Get the Feature type from df
                source=self._sources.get("columns"),
                exclude=self._exclude.get("columns"),
                **self._kwargs,  # type: ignore
            )

    def add_new_from(self, key: str, organism: str | None = None, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        self._kwargs.update({"organism": organism} if organism else {})
        self._update_registry(key, validated_only=False, **self._kwargs, **kwargs)

    def _add_new_from_columns(self, organism: str | None = None, **kwargs):
        """Deprecated to run by default during init."""
        warnings.warn(
            "`.add_new_from_columns()` is deprecated and will be removed in a future version. It's run by default during initialization.",
            DeprecationWarning,
            stacklevel=2,
        )
        pass

    def _replace_synonyms(
        self, key: str, syn_mapper: dict, values: pd.Series | pd.Index
    ):
        # replace the values in df
        std_values = values.map(lambda unstd_val: syn_mapper.get(unstd_val, unstd_val))
        # remove the standardized values from self.non_validated
        non_validated = [i for i in self.non_validated[key] if i not in syn_mapper]
        if len(non_validated) == 0:
            self._non_validated.pop(key, None)  # type: ignore
        else:
            self._non_validated[key] = non_validated  # type: ignore
        # logging
        n = len(syn_mapper)
        if n > 0:
            syn_mapper_print = _format_values(
                [f'"{k}" → "{v}"' for k, v in syn_mapper.items()], sep=""
            )
            s = "s" if n > 1 else ""
            logger.success(
                f'standardized {n} synonym{s} in "{key}": {colors.green(syn_mapper_print)}'
            )
        return std_values

    def standardize(self, key: str) -> None:
        """Replace synonyms with standardized values.

        Modifies the input dataset inplace.

        Args:
            key: The key referencing the column in the DataFrame to standardize.
        """
        # list is needed to avoid RuntimeError: dictionary changed size during iteration
        avail_keys = list(self.non_validated.keys())
        if len(avail_keys) == 0:
            logger.warning("values are already standardized")
            return

        if key == "all":
            for k in avail_keys:
                if k in self._fields:  # needed to exclude var_index
                    syn_mapper = standardize_categories(
                        self.non_validated[k],
                        field=self._fields[k],
                        using_key=self._using_key,
                        source=self._sources.get(k),
                        **self._kwargs,
                    )
                    self._df[k] = self._replace_synonyms(k, syn_mapper, self._df[k])
        else:
            if key not in avail_keys:
                if key in self._fields:
                    logger.info(f"No unstandardized values found for {key!r}")
                else:
                    raise KeyError(
                        f"{key!r} is not a valid key, available keys are: {_format_values(avail_keys)}!"
                    )
            else:
                if key in self._fields:  # needed to exclude var_index
                    syn_mapper = standardize_categories(
                        self.non_validated[key],
                        field=self._fields[key],
                        using_key=self._using_key,
                        source=self._sources.get(key),
                        **self._kwargs,
                    )
                    self._df[key] = self._replace_synonyms(
                        key, syn_mapper, self._df[key]
                    )

    def _update_registry(
        self, categorical: str, validated_only: bool = True, **kwargs
    ) -> None:
        if categorical == "all":
            self._update_registry_all(validated_only=validated_only, **kwargs)
        else:
            if categorical not in self.fields:
                raise ValidationError(
                    f"Feature {categorical} is not part of the fields!"
                )
            update_registry(
                values=_flatten_unique(self._df[categorical]),
                field=self.fields[categorical],
                key=categorical,
                using_key=self._using_key,
                validated_only=validated_only,
                source=self._sources.get(categorical),
                exclude=self._exclude.get(categorical),
                **kwargs,
            )
            # adding new records removes them from non_validated
            if not validated_only and self._non_validated:
                self._non_validated.pop(categorical, None)  # type: ignore

    def _update_registry_all(self, validated_only: bool = True, **kwargs):
        """Save labels for all features."""
        for name in self.fields.keys():
            self._update_registry(name, validated_only=validated_only, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate variables and categorical observations.

        This method also registers the validated records in the current instance:
        - from public sources
        - from the using_key instance

        Args:
            organism: The organism name.

        Returns:
            Whether the DataFrame is validated.
        """
        self._kwargs.update({"organism": organism} if organism else {})

        # add all validated records to the current instance
        self._update_registry_all()

        self._validated, self._non_validated = validate_categories_in_df(  # type: ignore
            self._df,
            fields=self.fields,
            using_key=self._using_key,
            sources=self._sources,
            exclude=self._exclude,
            **self._kwargs,
        )
        return self._validated

    def save_artifact(
        self,
        description: str | None = None,
        key: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Save the validated DataFrame and metadata.

        Args:
            description: Description of the DataFrame object.
            key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`.
                Artifacts with the same key form a revision family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        from lamindb.core._settings import settings

        if not self._validated:
            self.validate()
            if not self._validated:
                raise ValidationError("Dataset does not validate. Please curate.")

        # Make sure all labels are saved in the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            self._artifact = save_artifact(
                self._df,
                description=description,
                fields=self.fields,
                columns_field=self._columns_field,
                key=key,
                revises=revises,
                run=run,
                **self._kwargs,
            )
        finally:
            settings.verbosity = verbosity

        return self._artifact

    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't save any outputs."""
        from lamindb.core._context import context

        if context.run is not None:
            Run.filter(transform=context.run.transform, output_artifacts=None).exclude(
                uid=context.run.uid
            ).delete()


class AnnDataCurator(DataFrameCurator):
    """Curation flow for ``AnnData``.

    See also :class:`~lamindb.Curator`.

    Note that if genes are removed from the AnnData object, the object should be recreated using :meth:`~lamindb.Curator.from_anndata`.

    See :doc:`docs:cellxgene-curate` for instructions on how to curate against a specific cellxgene schema version.

    Args:
        data: The AnnData object or an AnnData-like path.
        var_index: The registry field for mapping the ``.var`` index.
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
        obs_columns: The registry field for mapping the ``.obs.columns``.
        using_key: A reference LaminDB instance.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping ``.obs.columns`` to Source records.
        exclude: A dictionary mapping column names to values to exclude from validation.
            When specific :class:`~bionty.Source` instances are pinned and may lack default values (e.g., "unknown" or "na"),
            using the exclude parameter ensures they are not validated.

    Examples:
        >>> import bionty as bt
        >>> curator = ln.Curator.from_anndata(
        ...     adata,
        ...     var_index=bt.Gene.ensembl_gene_id,
        ...     categoricals={
        ...         "cell_type_ontology_id": bt.CellType.ontology_id,
        ...         "donor_id": ln.ULabel.name
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        data: ad.AnnData | UPathStr,
        var_index: FieldAttr,
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict | None = None,
    ) -> None:
        from lamindb_setup.core import upath

        if isinstance(var_index, str):
            raise TypeError("var_index parameter has to be a bionty field")

        from ._artifact import data_is_anndata

        if sources is None:
            sources = {}
        if not data_is_anndata(data):
            raise TypeError(
                "data has to be an AnnData object or a path to AnnData-like"
            )
        if isinstance(data, ad.AnnData):
            self._adata = data
        else:  # pragma: no cover
            from lamindb.core.storage._backed_access import backed_access

            self._adata = backed_access(upath.create_path(data))

        if "symbol" in str(var_index):
            logger.warning(
                "indexing datasets with gene symbols can be problematic: https://docs.lamin.ai/faq/symbol-mapping"
            )

        self._data = data
        self._var_field = var_index
        super().__init__(
            df=self._adata.obs,
            categoricals=categoricals,
            columns=obs_columns,
            using_key=using_key,
            verbosity=verbosity,
            organism=organism,
            sources=sources,
            exclude=exclude,
            check_valid_keys=False,
        )
        self._obs_fields = categoricals or {}
        self._check_valid_keys(extra={"var_index"})

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(
        self, using_key: str | None = None, public: bool = False
    ) -> CurateLookup:
        """Lookup categories.

        Args:
            using_key: The instance where the lookup is performed.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, "var_index": self._var_field},
            using_key=using_key or self._using_key,
            public=public,
        )

    def _save_from_var_index(
        self, validated_only: bool = True, organism: str | None = None
    ):
        """Save variable records."""
        update_registry(
            values=list(self._adata.var.index),
            field=self.var_index,
            key="var_index",
            using_key=self._using_key,
            validated_only=validated_only,
            organism=organism,
            source=self._sources.get("var_index"),
            exclude=self._exclude.get("var_index"),
        )

    def _update_registry_all(self, validated_only: bool = True, **kwargs):
        """Save labels for all features."""
        self._save_from_var_index(validated_only=validated_only, **self._kwargs)
        for name in self._obs_fields.keys():
            self._update_registry(name, validated_only=validated_only, **self._kwargs)

    def add_new_from_var_index(self, organism: str | None = None, **kwargs):
        """Update variable records.

        Args:
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._save_from_var_index(validated_only=False, **self._kwargs, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate categories.

        This method also registers the validated records in the current instance.

        Args:
            organism: The organism name.

        Returns:
            Whether the AnnData object is validated.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        if self._using_key is not None and self._using_key != "default":
            logger.important(
                f"validating metadata using registries of instance {colors.italic(self._using_key)}"
            )

        # add all validated records to the current instance
        self._update_registry_all()

        validated_var, non_validated_var = validate_categories(
            self._adata.var.index,
            field=self._var_field,
            key="var_index",
            using_key=self._using_key,
            source=self._sources.get("var_index"),
            hint_print=".add_new_from_var_index()",
            exclude=self._exclude.get("var_index"),
            **self._kwargs,  # type: ignore
        )
        validated_obs, non_validated_obs = validate_categories_in_df(
            self._adata.obs,
            fields=self.categoricals,
            using_key=self._using_key,
            sources=self._sources,
            exclude=self._exclude,
            **self._kwargs,
        )
        self._non_validated = non_validated_obs  # type: ignore
        if len(non_validated_var) > 0:
            self._non_validated["var_index"] = non_validated_var  # type: ignore
        self._validated = validated_var and validated_obs
        return self._validated

    def standardize(self, key: str):
        """Replace synonyms with standardized values.

        Args:
            key: The key referencing the slot in `adata.obs` from which to draw terms. Same as the key in `categoricals`.

                - If "var_index", standardize the var.index.
                - If "all", standardize all obs columns and var.index.

        Inplace modification of the dataset.
        """
        if key in self._adata.obs.columns or key == "all":
            # standardize obs columns
            super().standardize(key)
        # in addition to the obs columns, standardize the var.index
        if key == "var_index" or key == "all":
            syn_mapper = standardize_categories(
                self._adata.var.index,
                field=self.var_index,
                using_key=self._using_key,
                source=self._sources.get("var_index"),
                **self._kwargs,
            )
            if "var_index" in self._non_validated:  # type: ignore
                self._adata.var.index = self._replace_synonyms(
                    "var_index", syn_mapper, self._adata.var.index
                )

    def save_artifact(
        self,
        description: str | None = None,
        key: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Save the validated ``AnnData`` and metadata.

        Args:
            description: A description of the ``AnnData`` object.
            key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`.
                Artifacts with the same key form a revision family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        from lamindb.core._settings import settings

        if not self._validated:
            self.validate()
            if not self._validated:
                raise ValidationError("Dataset does not validate. Please curate.")
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            self._artifact = save_artifact(
                self._data,
                adata=self._adata,
                description=description,
                columns_field=self.var_index,
                fields=self.categoricals,
                key=key,
                revises=revises,
                run=run,
                **self._kwargs,
            )
        finally:
            settings.verbosity = verbosity
        return self._artifact


class MuDataCurator:
    """Curation flow for a ``MuData`` object.

    See also :class:`~lamindb.Curator`.

    Note that if genes or other measurements are removed from the MuData object,
    the object should be recreated using :meth:`~lamindb.Curator.from_mudata`.

    Args:
        mdata: The MuData object to curate.
        var_index: The registry field for mapping the ``.var`` index for each modality.
            For example:
            ``{"modality_1": bt.Gene.ensembl_gene_id, "modality_2": ln.CellMarker.name}``
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
            Use modality keys to specify categoricals for MuData slots such as `"rna:cell_type": bt.CellType.name"`.
        using_key: A reference LaminDB instance.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping ``.obs.columns`` to Source records.
        exclude: A dictionary mapping column names to values to exclude from validation.
            When specific :class:`~bionty.Source` instances are pinned and may lack default values (e.g., "unknown" or "na"),
            using the exclude parameter ensures they are not validated.

    Examples:
        >>> import bionty as bt
        >>> curator = ln.Curator.from_mudata(
        ...     mdata,
        ...     var_index={
        ...         "rna": bt.Gene.ensembl_gene_id,
        ...         "adt": ln.CellMarker.name
        ...     },
        ...     categoricals={
        ...         "cell_type_ontology_id": bt.CellType.ontology_id,
        ...         "donor_id": ln.ULabel.name
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        mdata: MuData,
        var_index: dict[str, FieldAttr],
        categoricals: dict[str, FieldAttr] | None = None,
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict | None = None,  # {modality: {field: [values]}}
    ) -> None:
        if sources is None:
            sources = {}
        self._sources = sources
        if exclude is None:
            exclude = {}
        self._exclude = exclude
        self._mdata = mdata
        self._kwargs = {"organism": organism} if organism else {}
        self._var_fields = var_index
        self._verify_modality(self._var_fields.keys())
        self._obs_fields = self._parse_categoricals(categoricals)
        self._modalities = set(self._var_fields.keys()) | set(self._obs_fields.keys())
        self._using_key = using_key
        self._verbosity = verbosity
        self._obs_df_curator = None
        if "obs" in self._modalities:
            self._obs_df_curator = DataFrameCurator(
                df=mdata.obs,
                columns=Feature.name,
                categoricals=self._obs_fields.get("obs", {}),
                using_key=using_key,
                verbosity=verbosity,
                sources=self._sources.get("obs"),
                exclude=self._exclude.get("obs"),
                check_valid_keys=False,
                **self._kwargs,
            )
        self._mod_adata_curators = {
            modality: AnnDataCurator(
                data=mdata[modality],
                var_index=var_index.get(modality),
                categoricals=self._obs_fields.get(modality),
                using_key=using_key,
                verbosity=verbosity,
                sources=self._sources.get(modality),
                exclude=self._exclude.get(modality),
                **self._kwargs,
            )
            for modality in self._modalities
            if modality != "obs"
        }
        self._non_validated = None

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_fields

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    @property
    def non_validated(self) -> dict[str, dict[str, list[str]]]:
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return self._non_validated

    def _verify_modality(self, modalities: Iterable[str]):
        """Verify the modality exists."""
        for modality in modalities:
            if modality not in self._mdata.mod.keys():
                raise ValidationError(f"modality '{modality}' does not exist!")

    def _parse_categoricals(self, categoricals: dict[str, FieldAttr]) -> dict:
        """Parse the categorical fields."""
        prefixes = {f"{k}:" for k in self._mdata.mod.keys()}
        obs_fields: dict[str, dict[str, FieldAttr]] = {}
        for k, v in categoricals.items():
            if k not in self._mdata.obs.columns:
                raise ValidationError(f"column '{k}' does not exist in mdata.obs!")
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

    def lookup(
        self, using_key: str | None = None, public: bool = False
    ) -> CurateLookup:
        """Lookup categories.

        Args:
            using_key: The instance where the lookup is performed.
                if "public", the lookup is performed on the public reference.
        """
        obs_fields = {}
        for mod, fields in self._obs_fields.items():
            for k, v in fields.items():
                if k == "obs":
                    obs_fields[k] = v
                else:
                    obs_fields[f"{mod}:{k}"] = v
        return CurateLookup(
            categoricals=obs_fields,
            slots={
                **{f"{k}_var_index": v for k, v in self._var_fields.items()},
            },
            using_key=using_key or self._using_key,
            public=public,
        )

    def add_new_from_columns(
        self,
        modality: str,
        column_names: list[str] | None = None,
        organism: str | None = None,
        **kwargs,
    ):
        """Update columns records."""
        warnings.warn(
            "`.add_new_from_columns()` is deprecated and will be removed in a future version. It's run by default during initialization.",
            DeprecationWarning,
            stacklevel=2,
        )

    def add_new_from_var_index(
        self, modality: str, organism: str | None = None, **kwargs
    ):
        """Update variable records.

        Args:
            modality: The modality name.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self._kwargs.update({"organism": organism} if organism else {})
        self._mod_adata_curators[modality].add_new_from_var_index(
            **self._kwargs, **kwargs
        )

    def _update_registry_all(self):
        """Update all registries."""
        if self._obs_df_curator is not None:
            self._obs_df_curator._update_registry_all(
                validated_only=True, **self._kwargs
            )
        for _, adata_curator in self._mod_adata_curators.items():
            adata_curator._update_registry_all(validated_only=True, **self._kwargs)

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
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        self._kwargs.update({"organism": organism} if organism else {})
        modality = modality or "obs"
        if modality in self._mod_adata_curators:
            adata_curator = self._mod_adata_curators[modality]
            adata_curator.add_new_from(key=key, **self._kwargs, **kwargs)
        if modality == "obs":
            self._obs_df_curator.add_new_from(key=key, **self._kwargs, **kwargs)

    def validate(self, organism: str | None = None) -> bool:
        """Validate categories."""
        from lamindb.core._settings import settings

        self._kwargs.update({"organism": organism} if organism else {})
        if self._using_key is not None and self._using_key != "default":
            logger.important(
                f"validating using registries of instance {colors.italic(self._using_key)}"
            )

        # add all validated records to the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "error"
            self._update_registry_all()
        finally:
            settings.verbosity = verbosity

        self._non_validated = {}  # type: ignore

        obs_validated = True
        if "obs" in self._modalities:
            logger.info('validating categoricals in "obs"...')
            obs_validated &= self._obs_df_curator.validate(**self._kwargs)
            self._non_validated["obs"] = self._obs_df_curator.non_validated  # type: ignore
            logger.print("")

        mods_validated = True
        for modality, adata_curator in self._mod_adata_curators.items():
            logger.info(f'validating categoricals in modality "{modality}"...')
            mods_validated &= adata_curator.validate(**self._kwargs)
            if len(adata_curator.non_validated) > 0:
                self._non_validated[modality] = adata_curator.non_validated  # type: ignore
            logger.print("")

        self._validated = obs_validated & mods_validated
        return self._validated

    def standardize(self, key: str, modality: str | None = None):
        """Replace synonyms with standardized values.

        Args:
            key: The key referencing the slot in the `MuData`.
            modality: The modality name.

        Inplace modification of the dataset.
        """
        modality = modality or "obs"
        if modality in self._mod_adata_curators:
            adata_curator = self._mod_adata_curators[modality]
            adata_curator.standardize(key=key)
        if modality == "obs":
            self._obs_df_curator.standardize(key=key)

    def save_artifact(
        self,
        description: str | None = None,
        key: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Save the validated ``MuData`` and metadata.

        Args:
            description: A description of the ``MuData`` object.
            key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        from lamindb.core._settings import settings

        if not self._validated:
            self.validate()
            if not self._validated:
                raise ValidationError("Dataset does not validate. Please curate.")
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            self._artifact = save_artifact(
                self._mdata,
                description=description,
                columns_field=self.var_index,
                fields=self.categoricals,
                key=key,
                revises=revises,
                run=run,
                **self._kwargs,
            )
        finally:
            settings.verbosity = verbosity
        return self._artifact


def _maybe_curation_keys_not_present(nonval_keys: list[str], name: str):
    if (n := len(nonval_keys)) > 0:
        s = "s" if n > 1 else ""
        are = "are" if n > 1 else "is"
        raise ValidationError(
            f"key{s} passed to {name} {are} not present: {colors.yellow(_format_values(nonval_keys))}"
        )


class SOMACurator(BaseCurator):
    """Curation flow for ``tiledbsoma``.

    See also :class:`~lamindb.Curator`.

    Args:
        experiment_uri: A local or cloud path to a `tiledbsoma.Experiment`.
        var_index: The registry fields for mapping the `.var` indices for measurements.
            Should be in the form `{"measurement name": ("var column", field)}`.
            These keys should be used in the flattened form (`'{measurement name}__{column name in .var}'`)
            in `.standardize` or `.add_new_from`, see the output of `.var_index`.
        categoricals: A dictionary mapping categorical `.obs` columns to a registry field.
        obs_columns: The registry field for mapping the names of the `.obs` columns.
        organism: The organism name.
        sources: A dictionary mapping `.obs` columns to Source records.
        exclude: A dictionary mapping column names to values to exclude from validation.
            When specific :class:`~bionty.Source` instances are pinned and may lack default values (e.g., "unknown" or "na"),
            using the exclude parameter ensures they are not validated.

    Examples:
        >>> import bionty as bt
        >>> curator = ln.Curator.from_tiledbsoma(
        ...     "./my_array_store.tiledbsoma",
        ...     var_index={"RNA": ("var_id", bt.Gene.symbol)},
        ...     categoricals={
        ...         "cell_type_ontology_id": bt.CellType.ontology_id,
        ...         "donor_id": ln.ULabel.name
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        experiment_uri: UPathStr,
        var_index: dict[str, tuple[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict[str, str | list[str]] | None = None,
        using_key: str | None = None,
    ):
        self._obs_fields = categoricals or {}
        self._var_fields = var_index
        self._columns_field = obs_columns
        self._experiment_uri = UPath(experiment_uri)
        self._organism = organism
        self._using_key = using_key
        self._sources = sources or {}
        self._exclude = exclude or {}

        self._validated: bool | None = False
        self._non_validated_values: dict[str, list] | None = None
        self._validated_values: dict[str, list] = {}
        # filled by _check_save_keys
        self._n_obs: int | None = None
        self._valid_obs_keys: list[str] | None = None
        self._valid_var_keys: list[str] | None = None
        self._var_fields_flat: dict[str, FieldAttr] | None = None
        self._check_save_keys()

    # check that the provided keys in var_index and categoricals are available in the store
    # and save features
    def _check_save_keys(self):
        from lamindb.core.storage._tiledbsoma import _open_tiledbsoma

        with _open_tiledbsoma(self._experiment_uri, mode="r") as experiment:
            experiment_obs = experiment.obs
            self._n_obs = len(experiment_obs)
            valid_obs_keys = [k for k in experiment_obs.keys() if k != "soma_joinid"]
            self._valid_obs_keys = valid_obs_keys

            valid_var_keys = []
            ms_list = []
            for ms in experiment.ms.keys():
                ms_list.append(ms)
                var_ms = experiment.ms[ms].var
                valid_var_keys += [
                    f"{ms}__{k}" for k in var_ms.keys() if k != "soma_joinid"
                ]
            self._valid_var_keys = valid_var_keys

        # check validity of keys in categoricals
        nonval_keys = []
        for obs_key in self._obs_fields.keys():
            if obs_key not in valid_obs_keys:
                nonval_keys.append(obs_key)
        _maybe_curation_keys_not_present(nonval_keys, "categoricals")

        # check validity of keys in var_index
        self._var_fields_flat = {}
        nonval_keys = []
        for ms_key in self._var_fields.keys():
            var_key, var_field = self._var_fields[ms_key]
            var_key_flat = f"{ms_key}__{var_key}"
            if var_key_flat not in valid_var_keys:
                nonval_keys.append(f"({ms_key}, {var_key})")
            else:
                self._var_fields_flat[var_key_flat] = var_field
        _maybe_curation_keys_not_present(nonval_keys, "var_index")

        # check validity of keys in sources and exclude
        valid_arg_keys = valid_obs_keys + valid_var_keys + ["columns"]
        for name, dct in (("sources", self._sources), ("exclude", self._exclude)):
            nonval_keys = []
            for arg_key in dct.keys():
                if arg_key not in valid_arg_keys:
                    nonval_keys.append(arg_key)
            _maybe_curation_keys_not_present(nonval_keys, name)

        # register obs columns' names
        register_columns = list(self._obs_fields.keys())
        organism = check_registry_organism(
            self._columns_field.field.model, self._organism
        ).get("organism")
        update_registry(
            values=register_columns,
            field=self._columns_field,
            key="columns",
            using_key=self._using_key,
            validated_only=False,
            organism=organism,
            source=self._sources.get("columns"),
            exclude=self._exclude.get("columns"),
        )
        additional_columns = [k for k in valid_obs_keys if k not in register_columns]
        # no need to register with validated_only=True if columns are features
        if (
            len(additional_columns) > 0
            and self._columns_field.field.model is not Feature
        ):
            update_registry(
                values=additional_columns,
                field=self._columns_field,
                key="columns",
                using_key=self._using_key,
                validated_only=True,
                organism=organism,
                source=self._sources.get("columns"),
                exclude=self._exclude.get("columns"),
            )

    def validate(self):
        """Validate categories."""
        from lamindb.core.storage._tiledbsoma import _open_tiledbsoma

        validated = True
        self._non_validated_values = {}
        with _open_tiledbsoma(self._experiment_uri, mode="r") as experiment:
            for ms, (key, field) in self._var_fields.items():
                var_ms = experiment.ms[ms].var
                var_ms_key = f"{ms}__{key}"
                # it was already validated and cached
                if var_ms_key in self._validated_values:
                    continue
                var_ms_values = (
                    var_ms.read(column_names=[key]).concat()[key].to_pylist()
                )
                organism = check_registry_organism(
                    field.field.model, self._organism
                ).get("organism")
                update_registry(
                    values=var_ms_values,
                    field=field,
                    key=var_ms_key,
                    using_key=self._using_key,
                    validated_only=True,
                    organism=organism,
                    source=self._sources.get(var_ms_key),
                    exclude=self._exclude.get(var_ms_key),
                )
                _, non_val = validate_categories(
                    values=var_ms_values,
                    field=field,
                    key=var_ms_key,
                    using_key=self._using_key,
                    organism=organism,
                    source=self._sources.get(var_ms_key),
                    exclude=self._exclude.get(var_ms_key),
                )
                if len(non_val) > 0:
                    validated = False
                    self._non_validated_values[var_ms_key] = non_val
                else:
                    self._validated_values[var_ms_key] = var_ms_values

            obs = experiment.obs
            for key, field in self._obs_fields.items():
                # already validated and cached
                if key in self._validated_values:
                    continue
                values = pa.compute.unique(
                    obs.read(column_names=[key]).concat()[key]
                ).to_pylist()
                organism = check_registry_organism(
                    field.field.model, self._organism
                ).get("organism")
                update_registry(
                    values=values,
                    field=field,
                    key=key,
                    using_key=self._using_key,
                    validated_only=True,
                    organism=organism,
                    source=self._sources.get(key),
                    exclude=self._exclude.get(key),
                )
                _, non_val = validate_categories(
                    values=values,
                    field=field,
                    key=key,
                    using_key=self._using_key,
                    organism=organism,
                    source=self._sources.get(key),
                    exclude=self._exclude.get(key),
                )
                if len(non_val) > 0:
                    validated = False
                    self._non_validated_values[key] = non_val
                else:
                    self._validated_values[key] = values
        self._validated = validated
        return self._validated

    def _non_validated_values_field(self, key: str) -> tuple[list, FieldAttr]:
        assert self._non_validated_values is not None  # noqa: S101

        if key in self._valid_obs_keys:
            field = self._obs_fields[key]
        elif key in self._valid_var_keys:
            ms = key.partition("__")[0]
            field = self._var_fields[ms][1]
        else:
            raise KeyError(f"key {key} is invalid!")
        values = self._non_validated_values.get(key, [])
        return values, field

    def add_new_from(self, key: str) -> None:
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the `tiledbsoma` store.
                It should be `'{measurement name}__{column name in .var}'` for columns in `.var`
                or a column name in `.obs`.
        """
        if self._non_validated_values is None:
            raise ValidationError("Run .validate() first.")
        if key == "all":
            keys = list(self._non_validated_values.keys())
        else:
            avail_keys = list(
                chain(self._non_validated_values.keys(), self._validated_values.keys())
            )
            if key not in avail_keys:
                raise KeyError(
                    f"'{key!r}' is not a valid key, available keys are: {_format_values(avail_keys + ['all'])}!"
                )
            keys = [key]
        for k in keys:
            values, field = self._non_validated_values_field(k)
            if len(values) == 0:
                continue
            organism = check_registry_organism(field.field.model, self._organism).get(
                "organism"
            )
            update_registry(
                values=values,
                field=field,
                key=k,
                using_key=self._using_key,
                validated_only=False,
                organism=organism,
                source=self._sources.get(k),
                exclude=self._exclude.get(k),
            )
            # update non-validated values list but keep the key there
            # it will be removed by .validate()
            if k in self._non_validated_values:
                self._non_validated_values[k] = []

    @property
    def non_validated(self) -> dict[str, list]:
        """Return the non-validated features and labels."""
        non_val = {k: v for k, v in self._non_validated_values.items() if v != []}
        return non_val

    @property
    def var_index(self) -> dict[str, FieldAttr]:
        """Return the registry fields with flattened keys to validate variables indices against."""
        return self._var_fields_flat

    @property
    def categoricals(self) -> dict[str, FieldAttr]:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(
        self, using_key: str | None = None, public: bool = False
    ) -> CurateLookup:
        """Lookup categories.

        Args:
            using_key: The instance where the lookup is performed.
                if "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, **self._var_fields_flat},
            using_key=using_key or self._using_key,
            public=public,
        )

    def standardize(self, key: str):
        """Replace synonyms with standardized values.

        Modifies the dataset inplace.

        Args:
            key: The key referencing the slot in the `tiledbsoma` store.
                It should be `'{measurement name}__{column name in .var}'` for columns in `.var`
                or a column name in `.obs`.
        """
        if len(self.non_validated) == 0:
            logger.warning("values are already standardized")
            return
        avail_keys = list(self._non_validated_values.keys())
        if key == "all":
            keys = avail_keys
        else:
            if key not in avail_keys:
                raise KeyError(
                    f"'{key!r}' is not a valid key, available keys are: {_format_values(avail_keys + ['all'])}!"
                )
            keys = [key]

        for k in keys:
            values, field = self._non_validated_values_field(k)
            if len(values) == 0:
                continue
            if k in self._valid_var_keys:
                ms, _, slot_key = k.partition("__")
                slot = lambda experiment: experiment.ms[ms].var  # noqa: B023
            else:
                slot = lambda experiment: experiment.obs
                slot_key = k
            # errors if public ontology and the model has no organism
            # has to be fixed in bionty
            organism = check_registry_organism(field.field.model, self._organism).get(
                "organism"
            )
            syn_mapper = standardize_categories(
                values=values,
                field=field,
                using_key=self._using_key,
                source=self._sources.get(k),
                organism=organism,
            )
            if (n_syn_mapper := len(syn_mapper)) == 0:
                continue

            from lamindb.core.storage._tiledbsoma import _open_tiledbsoma

            with _open_tiledbsoma(self._experiment_uri, mode="r") as experiment:
                value_filter = f"{slot_key} in {list(syn_mapper.keys())}"
                table = slot(experiment).read(value_filter=value_filter).concat()

            if len(table) == 0:
                continue

            df = table.to_pandas()
            # map values
            df[slot_key] = df[slot_key].map(
                lambda val: syn_mapper.get(val, val)  # noqa
            )
            # write the mapped values
            with _open_tiledbsoma(self._experiment_uri, mode="w") as experiment:
                slot(experiment).write(pa.Table.from_pandas(df, schema=table.schema))
            # update non_validated dict
            non_val_k = [
                nv for nv in self._non_validated_values[k] if nv not in syn_mapper
            ]
            self._non_validated_values[k] = non_val_k

            syn_mapper_print = _format_values(
                [f'"{m_k}" → "{m_v}"' for m_k, m_v in syn_mapper.items()], sep=""
            )
            s = "s" if n_syn_mapper > 1 else ""
            logger.success(
                f'standardized {n_syn_mapper} synonym{s} in "{k}": {colors.green(syn_mapper_print)}'
            )

    def save_artifact(
        self,
        description: str | None = None,
        key: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """Save the validated `tiledbsoma` store and metadata.

        Args:
            description: A description of the ``tiledbsoma`` store.
            key: A path-like key to reference artifact in default storage,
                e.g., `"myfolder/mystore.tiledbsoma"`. Artifacts with the same key form a revision family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        from lamindb.core._data import add_labels

        if not self._validated:
            self.validate()
            if not self._validated:
                raise ValidationError("Dataset does not validate. Please curate.")

        artifact = Artifact(
            self._experiment_uri,
            description=description,
            key=key,
            revises=revises,
            run=run,
        )
        artifact.n_observations = self._n_obs
        artifact._accessor = "tiledbsoma"
        artifact.save()

        feature_sets = {}
        if len(self._obs_fields) > 0:
            organism = check_registry_organism(
                self._columns_field.field.model, self._organism
            ).get("organism")
            feature_sets["obs"] = FeatureSet.from_values(
                values=list(self._obs_fields.keys()),
                field=self._columns_field,
                organism=organism,
                raise_validation_error=False,
            )
        for ms in self._var_fields:
            var_key, var_field = self._var_fields[ms]
            organism = check_registry_organism(
                var_field.field.model, self._organism
            ).get("organism")
            feature_sets[f"{ms}__var"] = FeatureSet.from_values(
                values=self._validated_values[f"{ms}__{var_key}"],
                field=var_field,
                organism=organism,
                raise_validation_error=False,
            )
        artifact._feature_sets = feature_sets

        feature_ref_is_name = _ref_is_name(self._columns_field)
        features = Feature.lookup().dict()
        for key, field in self._obs_fields.items():
            feature = features.get(key)
            registry = field.field.model
            organism = check_registry_organism(field.field.model, self._organism).get(
                "organism"
            )
            labels = registry.from_values(
                values=self._validated_values[key], field=field, organism=organism
            )
            if len(labels) == 0:
                continue
            if hasattr(registry, "_name_field"):
                label_ref_is_name = field.field.name == registry._name_field
                add_labels(
                    artifact,
                    records=labels,
                    feature=feature,
                    feature_ref_is_name=feature_ref_is_name,
                    label_ref_is_name=label_ref_is_name,
                    from_curator=True,
                )

        return artifact.save()


class Curator(BaseCurator):
    """Dataset curator.

    A `Curator` object makes it easy to save validated & annotated artifacts.

    Example:

    >>> curator = ln.Curator.from_df(
    >>>     df,
    >>>     # define validation criteria as mappings
    >>>     columns=ln.Feature.name,  # map column names
    >>>     categoricals={"perturbation": ln.ULabel.name},  # map categories
    >>> )
    >>> curator.validate()  # validate the data in df
    >>> artifact = curator.save_artifact(description="my RNA-seq")
    >>> artifact.describe()  # see annotations

    `curator.validate()` maps values within `df` according to the mapping criteria and logs validated & problematic values.

    If you find non-validated values, you have several options:

    - new values found in the data can be registered using :meth:`~lamindb.core.DataFrameCurator.add_new_from`
    - non-validated values can be accessed using :meth:`~lamindb.core.DataFrameCurator.non_validated` and addressed manually
    """

    @classmethod
    @doc_args(DataFrameCurator.__doc__)
    def from_df(
        cls,
        df: pd.DataFrame,
        categoricals: dict[str, FieldAttr] | None = None,
        columns: FieldAttr = Feature.name,
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
    ) -> DataFrameCurator:
        """{}"""  # noqa: D415
        return DataFrameCurator(
            df=df,
            categoricals=categoricals,
            columns=columns,
            using_key=using_key,
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
        obs_columns: FieldAttr = Feature.name,
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> AnnDataCurator:
        """{}"""  # noqa: D415
        return AnnDataCurator(
            data=data,
            var_index=var_index,
            categoricals=categoricals,
            obs_columns=obs_columns,
            using_key=using_key,
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
        using_key: str | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
    ) -> MuDataCurator:
        """{}"""  # noqa: D415
        return MuDataCurator(
            mdata=mdata,
            var_index=var_index,
            categoricals=categoricals,
            using_key=using_key,
            verbosity=verbosity,
            organism=organism,
        )

    @classmethod
    @doc_args(SOMACurator.__doc__)
    def from_tiledbsoma(
        cls,
        experiment_uri: UPathStr,
        var_index: dict[str, tuple[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        using_key: str | None = None,
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict[str, str | list[str]] | None = None,
    ) -> SOMACurator:
        """{}"""  # noqa: D415
        return SOMACurator(
            experiment_uri=experiment_uri,
            var_index=var_index,
            categoricals=categoricals,
            obs_columns=obs_columns,
            using_key=using_key,
            organism=organism,
            sources=sources,
            exclude=exclude,
        )


def get_registry_instance(registry: Record, using_key: str | None = None) -> Record:
    """Get a registry instance using a specific instance."""
    if using_key is not None and using_key != "default":
        return registry.using(using_key)
    return registry


def get_current_filter_kwargs(registry: type[Record], kwargs: dict) -> dict:
    """Make sure the source and organism are saved in the same database as the registry."""
    from lamindb.core._settings import settings

    db = registry.filter().db
    source = kwargs.get("source")
    organism = kwargs.get("organism")
    filter_kwargs = kwargs.copy()
    try:
        verbosity = settings.verbosity
        settings.verbosity = "error"
        if isinstance(organism, Record) and organism._state.db != "default":
            if db is None or db == "default":
                organism_default = copy.copy(organism)
                # save the organism record in the default database
                organism_default.save()
                filter_kwargs["organism"] = organism_default
        if isinstance(source, Record) and source._state.db != "default":
            if db is None or db == "default":
                source_default = copy.copy(source)
                # save the source record in the default database
                source_default.save()
                filter_kwargs["source"] = source_default
    finally:
        settings.verbosity = verbosity
    return filter_kwargs


def inspect_instance(
    values: Iterable[str],
    field: FieldAttr,
    registry: type[Record],
    exclude: str | list | None = None,
    **kwargs,
):
    """Inspect values using a registry."""
    # inspect exclude values in the default instance
    values = list(values)
    include_validated = []
    if exclude is not None:
        exclude = [exclude] if isinstance(exclude, str) else exclude
        exclude = [i for i in exclude if i in values]
        if len(exclude) > 0:
            # exclude values are validated without source and organism
            inspect_result_exclude = registry.inspect(exclude, field=field, mute=True)
            # if exclude values are validated, remove them from the values
            values = [i for i in values if i not in inspect_result_exclude.validated]
            include_validated = inspect_result_exclude.validated

    inspect_result = registry.inspect(values, field=field, mute=True, **kwargs)
    inspect_result._validated += include_validated
    inspect_result._non_validated = [
        i for i in inspect_result.non_validated if i not in include_validated
    ]

    return inspect_result


def check_registry_organism(registry: Record, organism: str | None = None) -> dict:
    """Check if a registry needs an organism and return the organism name."""
    if hasattr(registry, "organism_id"):
        import bionty as bt

        if organism is None and bt.settings.organism is None:
            raise ValidationError(
                f"{registry.__name__} registry requires an organism!\n"
                "      → please pass an organism name via organism="
            )
        return {"organism": organism or bt.settings.organism.name}
    return {}


def validate_categories(
    values: Iterable[str],
    field: FieldAttr,
    key: str,
    using_key: str | None = None,
    organism: str | None = None,
    source: Record | None = None,
    exclude: str | list | None = None,
    hint_print: str | None = None,
) -> tuple[bool, list]:
    """Validate ontology terms in a pandas series using LaminDB registries.

    Args:
        values: The values to validate.
        field: The field attribute.
        key: The key referencing the slot in the DataFrame.
        using_key: A reference LaminDB instance.
        organism: The organism name.
        source: The source record.
        exclude: Exclude specific values from validation.
        standardize: Whether to standardize the values.
        hint_print: The hint to print that suggests fixing non-validated values.
    """
    from lamindb._from_values import _format_values
    from lamindb.core._settings import settings

    model_field = f"{field.field.model.__name__}.{field.field.name}"

    def _log_mapping_info():
        logger.indent = ""
        logger.info(f'mapping "{key}" on {colors.italic(model_field)}')
        logger.indent = "  "

    registry = field.field.model

    # {"organism": organism_name/organism_record}
    kwargs = check_registry_organism(registry, organism)
    kwargs.update({"source": source} if source else {})
    kwargs_current = get_current_filter_kwargs(registry, kwargs)

    # inspect values from the default instance
    inspect_result = inspect_instance(
        values=values,
        field=field,
        registry=registry,
        exclude=exclude,
        **kwargs_current,
    )
    non_validated = inspect_result.non_validated
    syn_mapper = inspect_result.synonyms_mapper

    # inspect the non-validated values from the using_key instance
    values_validated = []
    if using_key is not None and using_key != "default" and non_validated:
        registry_using = get_registry_instance(registry, using_key)
        inspect_result = inspect_instance(
            values=non_validated,
            field=field,
            registry=registry_using,
            exclude=exclude,
            **kwargs,
        )
        non_validated = inspect_result.non_validated
        values_validated += inspect_result.validated
        syn_mapper.update(inspect_result.synonyms_mapper)

    # inspect the non-validated values from public (bionty only)
    if hasattr(registry, "public"):
        verbosity = settings.verbosity
        try:
            settings.verbosity = "error"
            public_records = registry.from_values(
                non_validated,
                field=field,
                **kwargs_current,
            )
            values_validated += [getattr(r, field.field.name) for r in public_records]
        finally:
            settings.verbosity = verbosity

    # logging messages
    non_validated_hint_print = hint_print or f'.add_new_from("{key}")'
    non_validated = [i for i in non_validated if i not in values_validated]
    n_non_validated = len(non_validated)
    if n_non_validated == 0:
        logger.indent = ""
        logger.success(f'"{key}" is validated against {colors.italic(model_field)}')
        return True, []
    else:
        are = "is" if n_non_validated == 1 else "are"
        s = "" if n_non_validated == 1 else "s"
        print_values = _format_values(non_validated)
        warning_message = f"{colors.red(f'{n_non_validated} term{s}')} {are} not validated: {colors.red(print_values)}\n"
        if syn_mapper:
            s = "" if len(syn_mapper) == 1 else "s"
            syn_mapper_print = _format_values(
                [f'"{k}" → "{v}"' for k, v in syn_mapper.items()], sep=""
            )
            hint_msg = f'.standardize("{key}")'
            warning_message += f"    {colors.yellow(f'{len(syn_mapper)} synonym{s}')} found: {colors.yellow(syn_mapper_print)}\n    → curate synonyms via {colors.cyan(hint_msg)}"
        if n_non_validated > len(syn_mapper):
            if syn_mapper:
                warning_message += "    for remaining terms:\n"
            warning_message += f"    → fix typos, remove non-existent values, or save terms via {colors.cyan(non_validated_hint_print)}"

        if logger.indent == "":
            _log_mapping_info()
        logger.warning(warning_message)
        logger.indent = ""
        return False, non_validated


def standardize_categories(
    values: Iterable[str],
    field: FieldAttr,
    using_key: str | None = None,
    organism: str | None = None,
    source: Record | None = None,
) -> dict:
    """Get a synonym mapper."""
    registry = field.field.model
    if not hasattr(registry, "standardize"):
        return {}
    # standardize values using the default instance
    syn_mapper = registry.standardize(
        values,
        field=field.field.name,
        organism=organism,
        source=source,
        mute=True,
        return_mapper=True,
    )

    if len(values) > len(syn_mapper):  # type: ignore
        # standardize values using the using_key instance
        if using_key is not None and using_key != "default":
            registry_using = get_registry_instance(registry, using_key)
            syn_mapper.update(
                registry_using.standardize(
                    [v for v in values if v not in syn_mapper],
                    field=field.field.name,
                    organism=organism,
                    source=source,
                    mute=True,
                    return_mapper=True,
                )
            )
    return syn_mapper


def validate_categories_in_df(
    df: pd.DataFrame,
    fields: dict[str, FieldAttr],
    using_key: str | None = None,
    sources: dict[str, Record] = None,
    exclude: dict | None = None,
    **kwargs,
) -> tuple[bool, dict]:
    """Validate categories in DataFrame columns using LaminDB registries."""
    if not fields:
        return True, {}

    if sources is None:
        sources = {}
    validated = True
    non_validated = {}
    for key, field in fields.items():
        is_val, non_val = validate_categories(
            df[key],
            field=field,
            key=key,
            using_key=using_key,
            source=sources.get(key),
            exclude=exclude.get(key) if exclude else None,
            **kwargs,
        )
        validated &= is_val
        if len(non_val) > 0:
            non_validated[key] = non_val
    return validated, non_validated


def save_artifact(
    data: pd.DataFrame | ad.AnnData | MuData,
    fields: dict[str, FieldAttr] | dict[str, dict[str, FieldAttr]],
    columns_field: FieldAttr | dict[str, FieldAttr],
    description: str | None = None,
    organism: str | None = None,
    adata: ad.AnnData | None = None,
    key: str | None = None,
    revises: Artifact | None = None,
    run: Run | None = None,
) -> Artifact:
    """Save all metadata with an Artifact.

    Args:
        data: The DataFrame or AnnData object to save.
        fields: A dictionary mapping obs_column to registry_field.
        columns_field: The registry field to validate variables index against.
        description: A description of the artifact.
        organism: The organism name.
        adata: The AnnData object to save and get n_observations, must be provided if data is a path.
        type: The artifact type.
        key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
        revises: Previous version of the artifact. Triggers a revision.
        run: The run that creates the artifact.

    Returns:
        The saved Artifact.
    """
    from ._artifact import data_is_anndata
    from .core._data import add_labels

    artifact = None
    if data_is_anndata(data):
        assert adata is not None  # noqa: S101
        artifact = Artifact.from_anndata(
            data, description=description, key=key, revises=revises, run=run
        )
        artifact.n_observations = adata.shape[0]
        data = adata

    elif isinstance(data, pd.DataFrame):
        artifact = Artifact.from_df(
            data, description=description, key=key, revises=revises, run=run
        )
    else:
        try:
            from mudata import MuData

            if isinstance(data, MuData):
                artifact = Artifact.from_mudata(
                    data,
                    description=description,
                    key=key,
                    revises=revises,
                    run=run,
                )
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

    def _add_labels(
        data,
        artifact: Artifact,
        fields: dict[str, FieldAttr],
        feature_ref_is_name: bool | None = None,
    ):
        features = Feature.lookup().dict()
        for key, field in fields.items():
            feature = features.get(key)
            registry = field.field.model
            filter_kwargs = check_registry_organism(registry, organism)
            filter_kwargs_current = get_current_filter_kwargs(registry, filter_kwargs)
            df = data if isinstance(data, pd.DataFrame) else data.obs
            labels = registry.from_values(
                df[key],
                field=field,
                **filter_kwargs_current,
            )
            if len(labels) == 0:
                continue
            label_ref_is_name = None
            if hasattr(registry, "_name_field"):
                label_ref_is_name = field.field.name == registry._name_field
            add_labels(
                artifact,
                records=labels,
                feature=feature,
                feature_ref_is_name=feature_ref_is_name,
                label_ref_is_name=label_ref_is_name,
                from_curator=True,
            )

    if artifact._accessor == "MuData":
        for modality, modality_fields in fields.items():
            column_field_modality = columns_field.get(modality)
            if modality == "obs":
                _add_labels(
                    data,
                    artifact,
                    modality_fields,
                    feature_ref_is_name=(
                        None
                        if column_field_modality is None
                        else _ref_is_name(column_field_modality)
                    ),
                )
            else:
                _add_labels(
                    data[modality],
                    artifact,
                    modality_fields,
                    feature_ref_is_name=(
                        None
                        if column_field_modality is None
                        else _ref_is_name(column_field_modality)
                    ),
                )
    else:
        _add_labels(
            data, artifact, fields, feature_ref_is_name=_ref_is_name(columns_field)
        )

    slug = ln_setup.settings.instance.slug
    if ln_setup.settings.instance.is_remote:  # pragma: no cover
        logger.important(f"go to https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact


def _flatten_unique(series: pd.Series[list[Any] | Any]) -> list[Any]:
    """Flatten a Pandas series containing lists or single items into a unique list of elements."""
    result = set()

    for item in series:
        if isinstance(item, list):
            result.update(item)
        else:
            result.add(item)

    return list(result)


def update_registry(
    values: list[str],
    field: FieldAttr,
    key: str,
    using_key: str | None = None,
    validated_only: bool = True,
    df: pd.DataFrame | None = None,
    organism: str | None = None,
    dtype: str | None = None,
    source: Record | None = None,
    exclude: str | list | None = None,
    **kwargs,
) -> None:
    """Save features or labels records in the default instance from the using_key instance.

    Args:
        values: A list of values to be saved as labels.
        field: The FieldAttr object representing the field for which labels are being saved.
        key: The name of the feature to save.
        using_key: The name of the instance from which to transfer labels (if applicable).
        validated_only: If True, only save validated labels.
        df: A DataFrame to save labels from.
        organism: The organism name.
        dtype: The type of the feature.
        source: The source record.
        exclude: Values to exclude from inspect.
        kwargs: Additional keyword arguments to pass to the registry model to create new records.
    """
    from lamindb._save import save as ln_save
    from lamindb.core._settings import settings

    registry = field.field.model
    filter_kwargs = check_registry_organism(registry, organism)
    filter_kwargs.update({"source": source} if source else {})
    if not values:
        return

    verbosity = settings.verbosity
    try:
        settings.verbosity = "error"
        labels_saved: dict = {"from public": [], "new": []}

        # inspect the default instance and save validated records from public
        filter_kwargs_current = get_current_filter_kwargs(registry, filter_kwargs)
        existing_and_public_records = registry.from_values(
            list(values), field=field, **filter_kwargs_current
        )
        existing_and_public_labels = [
            getattr(r, field.field.name) for r in existing_and_public_records
        ]
        # public records that are not already in the database
        public_records = [r for r in existing_and_public_records if r._state.adding]
        # here we check to only save the public records if they are from the specified source
        # we check the uid because r.source and source can be from different instances
        if source:
            public_records = [r for r in public_records if r.source.uid == source.uid]
        if len(public_records) > 0:
            settings.verbosity = "info"
            logger.info(f"saving validated records of '{key}'")
            settings.verbosity = "error"
            ln_save(public_records)
            labels_saved["from public"] = [
                getattr(r, field.field.name) for r in public_records
            ]
        # non-validated records from the default instance
        non_validated_labels = [
            i for i in values if i not in existing_and_public_labels
        ]

        # inspect and save validated records the using_key instance
        (
            labels_saved[f"from {using_key}"],
            non_validated_labels,
        ) = update_registry_from_using_instance(
            non_validated_labels,
            field=field,
            using_key=using_key,
            exclude=exclude,
            **filter_kwargs,
        )

        # save non-validated/new records
        labels_saved["new"] = non_validated_labels
        if not validated_only:
            non_validated_records = []
            if df is not None and registry == Feature:
                nonval_columns = Feature.inspect(df.columns, mute=True).non_validated
                non_validated_records = Feature.from_df(df.loc[:, nonval_columns])
            else:
                if "organism" in filter_kwargs:
                    # make sure organism record is saved to the current instance
                    filter_kwargs["organism"] = _save_organism(name=organism)
                init_kwargs = {}
                for value in labels_saved["new"]:
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

        # save parent labels for ulabels, for example a parent label "project" for label "project001"
        if registry == ULabel and field.field.name == "name":
            save_ulabels_parent(values, field=field, key=key)

    finally:
        settings.verbosity = verbosity

    log_saved_labels(
        labels_saved,
        key=key,
        model_field=f"{registry.__name__}.{field.field.name}",
        validated_only=validated_only,
    )


def log_saved_labels(
    labels_saved: dict,
    key: str,
    model_field: str,
    validated_only: bool = True,
) -> None:
    """Log the saved labels."""
    from ._from_values import _format_values

    model_field = colors.italic(model_field)
    for k, labels in labels_saved.items():
        if not labels:
            continue
        if k == "new" and validated_only:
            continue
        else:
            k = "" if k == "new" else f"{colors.green(k)} "
            # the term "transferred" stresses that this is always in the context of transferring
            # labels from a public ontology or a different instance to the present instance
            s = "s" if len(labels) > 1 else ""
            logger.success(
                f'added {len(labels)} record{s} {k}with {model_field} for "{key}": {_format_values(labels)}'
            )


def save_ulabels_parent(values: list[str], field: FieldAttr, key: str) -> None:
    """Save a parent label for the given labels."""
    registry = field.field.model
    assert registry == ULabel  # noqa: S101
    all_records = registry.from_values(list(values), field=field)
    is_feature = registry.filter(name=f"{key}").one_or_none()
    if is_feature is None:
        is_feature = registry(name=f"{key}").save()
        logger.important(f"Created a parent ULabel: {is_feature}")
    is_feature.children.add(*all_records)


def update_registry_from_using_instance(
    values: list[str],
    field: FieldAttr,
    using_key: str | None = None,
    exclude: str | list | None = None,
    **kwargs,
) -> tuple[list[str], list[str]]:
    """Save features or labels records from the using_key instance.

    Args:
        values: A list of values to be saved as labels.
        field: The FieldAttr object representing the field for which labels are being saved.
        using_key: The name of the instance from which to transfer labels (if applicable).
        kwargs: Additional keyword arguments to pass to the registry model.

    Returns:
        A tuple containing the list of saved labels and the list of non-saved labels.
    """
    labels_saved = []
    not_saved = values

    if using_key is not None and using_key != "default":
        registry_using = get_registry_instance(field.field.model, using_key)

        inspect_result_using = inspect_instance(
            values=values,
            field=field,
            registry=registry_using,
            exclude=exclude,
            **kwargs,
        )
        labels_using = registry_using.filter(
            **{f"{field.field.name}__in": inspect_result_using.validated}
        ).all()
        for label_using in labels_using:
            label_using.save()
            labels_saved.append(getattr(label_using, field.field.name))
        not_saved = inspect_result_using.non_validated

    return labels_saved, not_saved


def _save_organism(name: str):
    """Save an organism record."""
    import bionty as bt

    organism = bt.Organism.filter(name=name).one_or_none()
    if organism is None:
        organism = bt.Organism.from_source(name=name)
        if organism is None:
            raise ValidationError(
                f'Organism "{name}" not found\n'
                f'      → please save it: bt.Organism(name="{name}").save()'
            )
        organism.save()
    return organism


def _ref_is_name(field: FieldAttr) -> bool | None:
    """Check if the reference field is a name field."""
    from ._can_curate import get_name_field

    name_field = get_name_field(field.field.model)
    return field.field.name == name_field


Curate = Curator  # backward compat
