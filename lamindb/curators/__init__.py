"""Curators.

.. autosummary::
   :toctree: .

   DataFrameCurator

"""

from __future__ import annotations

import copy
import random
import re
from importlib import resources
from itertools import chain
from typing import TYPE_CHECKING, Any, Literal

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
import pandera as pda
import pyarrow as pa
from lamin_utils import colors, logger
from lamindb_setup.core import deprecated, upath
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.upath import UPath

from lamindb.core.storage._backed_access import backed_access

from ._cellxgene_schemas import _read_schema_versions

if TYPE_CHECKING:
    from anndata import AnnData
    from lamindb_setup.core.types import UPathStr

    from lamindb.base.types import FieldAttr
    from lamindb.models import Record
from lamindb._feature import parse_dtype, parse_dtype_single_cat
from lamindb.base.types import FieldAttr  # noqa
from lamindb.core._data import add_labels
from lamindb.core._feature_manager import parse_staged__schemas_m2m_from_anndata
from lamindb.core._settings import settings
from lamindb.models import (
    Artifact,
    CanCurate,
    Collection,
    Feature,
    Record,
    Run,
    Schema,
    ULabel,
)

from .._artifact import data_is_anndata
from .._from_values import _format_values
from ..errors import InvalidArgument, ValidationError

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableMapping
    from typing import Any

    from lamindb_setup.core.types import UPathStr
    from mudata import MuData
    from spatialdata import SpatialData

    from lamindb._query_set import RecordList


def strip_ansi_codes(text):
    # This pattern matches ANSI escape sequences
    ansi_pattern = re.compile(r"\x1b\[[0-9;]*m")
    return ansi_pattern.sub("", text)


class CurateLookup:
    """Lookup categories from the reference instance.

    Args:
        categoricals: A dictionary of categorical fields to lookup.
        slots: A dictionary of slot fields to lookup.
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
        public: bool = False,
    ) -> None:
        slots = slots or {}
        self._categoricals = {**categoricals, **slots}
        self._public = public

    def __getattr__(self, name):
        if name in self._categoricals:
            registry = self._categoricals[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public().lookup()
            else:
                return registry.lookup()
        raise AttributeError(
            f'"{self.__class__.__name__}" object has no attribute "{name}"'
        )

    def __getitem__(self, name):
        if name in self._categoricals:
            registry = self._categoricals[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public().lookup()
            else:
                return registry.lookup()
        raise AttributeError(
            f'"{self.__class__.__name__}" object has no attribute "{name}"'
        )

    def __repr__(self) -> str:
        if len(self._categoricals) > 0:
            getattr_keys = "\n ".join(
                [f".{key}" for key in self._categoricals if key.isidentifier()]
            )
            getitem_keys = "\n ".join(
                [str([key]) for key in self._categoricals if not key.isidentifier()]
            )
            ref = "public" if self._public else "registries"
            return (
                f"Lookup objects from the {colors.italic(ref)}:\n "
                f"{colors.green(getattr_keys)}\n "
                f"{colors.green(getitem_keys)}\n"
                'Example:\n    → categories = curator.lookup()["cell_type"]\n'
                "    → categories.alveolar_type_1_fibroblast_cell\n\n"
                "To look up public ontologies, use .lookup(public=True)"
            )
        else:  # pdagma: no cover
            return colors.warning("No fields are found!")


class Curator:
    """Dataset curator.

    A `Curator` object makes it easy to validate, standardize & annotate datasets.

    Example:

    >>> curator = ln.Curator(
    >>>     dataset,
    >>>     # define validation criteria as mappings
    >>>     columns=Feature.name,  # map column names
    >>>     categoricals={"perturbation": ULabel.name},  # map categories
    >>> )
    >>> curator.validate()  # validate the data in df
    >>> artifact = curator.save_artifact(description="my RNA-seq")
    >>> artifact.describe()  # see annotations

    `curator.validate()` maps values within `df` according to the mapping criteria and logs validated & problematic values.

    If you find non-validated values, you have several options:

    - new values found in the data can be registered using :meth:`~lamindb.core.DataFrameCatCurator.add_new_from`
    - non-validated values can be accessed using :meth:`~lamindb.core.DataFrameCatCurator.non_validated` and addressed manually
    """

    def __init__(self, dataset: Any, schema: Schema | None = None):
        self._dataset: Any = dataset  # pass the dataset as a UPathStr or data object
        self._artifact: Artifact = None  # pass the dataset as a non-curated artifact
        self._schema: Schema | None = schema
        self._is_validated: bool = False
        self._cat_curator: CatCurator = None  # is None for CatCurator curators

    def validate(self) -> bool | str:
        """Validate dataset.

        This method also registers the validated records in the current instance.

        Returns:
            The boolean `True` if the dataset is validated. Otherwise, a string with the error message.
        """
        pass  # pdagma: no cover

    def standardize(self, key: str) -> None:
        """Replace synonyms with standardized values.

        Inplace modification of the dataset.

        Args:
            key: The name of the column to standardize.

        Returns:
            None
        """
        pass  # pdagma: no cover

    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        # Note that this docstring has to be consistent with the Artifact()
        # constructor signature
        """Save an annotated artifact.

        Args:
            key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
            description: A description.
            revises: Previous version of the artifact. Is an alternative way to passing `key` to trigger a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        if not self._is_validated:
            self.validate()  # raises ValidationError if doesn't validate
        return self._cat_curator.save_artifact(
            key=key, description=description, revises=revises, run=run
        )


class DataFrameCurator(Curator):
    """Curator for a DataFrame object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    Args:
        dataset: The DataFrame-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Examples:
        >>> curator = curators.DataFrameCurator(
        ...     df,
        ...     schema,
        ... )
    """

    def __init__(
        self,
        dataset: pd.DataFrame,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if schema.n > 0:
            # populate features
            non_categoricals = {}
            categoricals = {}
            for feature in schema.features.all():
                pda_dtype = (
                    feature.dtype if not feature.dtype.startswith("cat") else "category"
                )
                non_categoricals[feature.name] = pda.Column(pda_dtype)
                if feature.dtype.startswith("cat"):
                    categoricals[feature.name] = parse_dtype(feature.dtype)[0]["field"]
            self._pda_schema = pda.DataFrameSchema(
                non_categoricals, coerce=schema.coerce_dtype
            )
            # now deal with categorical features using the old-style curator
            self._cat_curator = DataFrameCatCurator(
                dataset,
                categoricals=categoricals,
            )
        else:
            assert schema.itype is not None  # noqa: S101

    def validate(self) -> None:
        if self._schema.n > 0:
            self._cat_curator.validate()
            try:
                self._pda_schema.validate(self._dataset)
                if self._cat_curator._is_validated:
                    self._is_validated = True
                else:
                    self._is_validated = False
                    raise ValidationError(
                        self._cat_curator._validate_category_error_messages
                    )
            except pda.errors.SchemaError as err:
                self._is_validated = False
                # .exconly() doesn't exist on SchemaError
                raise ValidationError(str(err)) from err
        else:
            result = parse_dtype_single_cat(self._schema.itype)
            registry: CanCurate = result["registry"]
            inspector = registry.inspect(
                self._dataset.columns, result["field"], mute=True
            )
            if len(inspector.non_validated):
                self._is_validated = False
                raise ValidationError(
                    f"Invalid column identifiers found: {inspector.non_validated}"
                )


class AnnDataCurator(Curator):
    """Curator for a DataFrame object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    Args:
        dataset: The AnnData-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Examples:
        >>> curator = curators.DataFrameCurator(
        ...     df,
        ...     schema,
        ... )
    """

    def __init__(
        self,
        dataset: AnnData,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_anndata(dataset):
            raise InvalidArgument("dataset must be AnnData-like.")
        if schema.otype != "AnnData":
            raise InvalidArgument("Schema otype must be 'AnnData'.")
        self._obs_curator = DataFrameCurator(dataset.obs, schema._get_component("obs"))
        self._var_curator = DataFrameCurator(
            dataset.var.T, schema._get_component("var")
        )

    def validate(self) -> None:
        self._obs_curator.validate()
        self._var_curator.validate()
        self._is_validated = True

    def save_artifact(self, *, key=None, description=None, revises=None, run=None):
        result = parse_dtype_single_cat(self._var_curator._schema.itype)
        return save_artifact(  # type: ignore
            self._dataset,
            description=description,
            fields=self._obs_curator._cat_curator.categoricals,
            columns_field=result["field"],
            key=key,
            revises=revises,
            run=run,
            schema=self,
        )


class CatCurator(Curator):
    def __init__(
        self, *, dataset, categoricals, sources, organism, exclude, columns_field=None
    ):
        super().__init__(dataset=dataset)
        self._categoricals = categoricals or {}
        self._non_validated = None
        self._organism = organism
        self._sources = sources or {}
        self._exclude = exclude or {}
        self._columns_field = columns_field
        self._validate_category_error_messages: str = ""

    @property
    def non_validated(self) -> dict[str, list[str]]:
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return self._non_validated

    @property
    def categoricals(self) -> dict:
        """Return the columns fields to validate against."""
        return self._categoricals

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

    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        from lamindb.core._settings import settings

        if not self._is_validated:
            self.validate()  # returns True or False
            if not self._is_validated:  # need to raise error manually
                raise ValidationError("Dataset does not validate. Please curate.")

        # Make sure all labels are saved in the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"
            self._artifact = save_artifact(  # type: ignore
                self._dataset,
                description=description,
                fields=self.categoricals,
                columns_field=self._columns_field,
                key=key,
                revises=revises,
                run=run,
                schema=None,
                organism=self._organism,
            )
        finally:
            settings.verbosity = verbosity

        return self._artifact


class DataFrameCatCurator(CatCurator):
    """Curation flow for a DataFrame object.

    See also :class:`~lamindb.Curator`.

    Args:
        df: The DataFrame object to curate.
        columns: The field attribute for the feature column.
        categoricals: A dictionary mapping column names to registry_field.
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
        ...         "donor_id": ULabel.name
        ...     }
        ... )
    """

    def __init__(
        self,
        df: pd.DataFrame,
        columns: FieldAttr = Feature.name,
        categoricals: dict[str, FieldAttr] | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict | None = None,
    ) -> None:
        from lamindb.core._settings import settings

        if organism is not None and not isinstance(organism, str):
            raise ValueError("organism must be a string such as 'human' or 'mouse'!")

        settings.verbosity = verbosity
        self._non_validated = None
        super().__init__(
            dataset=df,
            columns_field=columns,
            organism=organism,
            categoricals=categoricals,
            sources=sources,
            exclude=exclude,
        )
        self._save_columns()

    def lookup(self, public: bool = False) -> CurateLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._categoricals,
            slots={"columns": self._columns_field},
            public=public,
        )

    def _save_columns(self, validated_only: bool = True) -> None:
        """Save column name records."""
        # Always save features specified as the fields keys
        update_registry(
            values=list(self.categoricals.keys()),
            field=self._columns_field,
            key="columns",
            validated_only=False,
            source=self._sources.get("columns"),
            exclude=self._exclude.get("columns"),
        )

        # Save the rest of the columns based on validated_only
        additional_columns = set(self._dataset.columns) - set(self.categoricals.keys())
        if additional_columns:
            update_registry(
                values=list(additional_columns),
                field=self._columns_field,
                key="columns",
                validated_only=validated_only,
                df=self._dataset,  # Get the Feature type from df
                source=self._sources.get("columns"),
                exclude=self._exclude.get("columns"),
            )

    @deprecated(new_name="is run by default")
    def add_new_from_columns(self, organism: str | None = None, **kwargs):
        pass

    def validate(self) -> bool:
        """Validate variables and categorical observations.

        This method also registers the validated records in the current instance:
        - from public sources

        Args:
            organism: The organism name.

        Returns:
            Whether the DataFrame is validated.
        """
        # add all validated records to the current instance
        self._update_registry_all()
        self._validate_category_error_messages = ""  # reset the error messages
        self._is_validated, self._non_validated = validate_categories_in_df(  # type: ignore
            self._dataset,
            fields=self.categoricals,
            sources=self._sources,
            exclude=self._exclude,
            curator=self,
            organism=self._organism,
        )
        return self._is_validated

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
                if k in self._categoricals:  # needed to exclude var_index
                    syn_mapper = standardize_categories(
                        self.non_validated[k],
                        field=self._categoricals[k],
                        source=self._sources.get(k),
                    )
                    self._dataset[k] = self._replace_synonyms(
                        k, syn_mapper, self._dataset[k]
                    )
        else:
            if key not in avail_keys:
                if key in self._categoricals:
                    logger.info(f"No unstandardized values found for {key!r}")
                else:
                    raise KeyError(
                        f"{key!r} is not a valid key, available keys are: {_format_values(avail_keys)}!"
                    )
            else:
                if key in self._categoricals:  # needed to exclude var_index
                    syn_mapper = standardize_categories(
                        self.non_validated[key],
                        field=self._categoricals[key],
                        source=self._sources.get(key),
                        organism=self._organism,
                    )
                    self._dataset[key] = self._replace_synonyms(
                        key, syn_mapper, self._dataset[key]
                    )

    def _update_registry_all(self, validated_only: bool = True, **kwargs):
        """Save labels for all features."""
        for name in self.categoricals.keys():
            self._update_registry(name, validated_only=validated_only, **kwargs)

    def _update_registry(
        self, categorical: str, validated_only: bool = True, **kwargs
    ) -> None:
        if categorical == "all":
            self._update_registry_all(validated_only=validated_only, **kwargs)
        else:
            if categorical not in self.categoricals:
                raise ValidationError(
                    f"Feature {categorical} is not part of the fields!"
                )
            update_registry(
                values=_flatten_unique(self._dataset[categorical]),
                field=self.categoricals[categorical],
                key=categorical,
                validated_only=validated_only,
                source=self._sources.get(categorical),
                exclude=self._exclude.get(categorical),
                organism=self._organism,
            )
            # adding new records removes them from non_validated
            if not validated_only and self._non_validated:
                self._non_validated.pop(categorical, None)  # type: ignore

    def add_new_from(self, key: str, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        self._update_registry(key, validated_only=False, **kwargs)

    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't save any outputs."""
        from lamindb.core._context import context

        if context.run is not None:
            Run.filter(transform=context.run.transform, output_artifacts=None).exclude(
                uid=context.run.uid
            ).delete()


class AnnDataCatCurator(CatCurator):
    """Manage categorical curation.

    Args:
        data: The AnnData object or an AnnData-like path.
        var_index: The registry field for mapping the ``.var`` index.
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
        obs_columns: The registry field for mapping the ``.obs.columns``.
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
        ...         "donor_id": ULabel.name
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
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict | None = None,
    ) -> None:
        from lamindb_setup.core import upath

        if isinstance(var_index, str):
            raise TypeError("var_index parameter has to be a bionty field")

        if sources is None:
            sources = {}
        if not data_is_anndata(data):
            raise TypeError(
                "data has to be an AnnData object or a path to AnnData-like"
            )
        if isinstance(data, ad.AnnData):
            self._adata = data
        else:  # pdagma: no cover
            from lamindb.core.storage._backed_access import backed_access

            self._adata = backed_access(upath.create_path(data))

        if "symbol" in str(var_index):
            logger.warning(
                "indexing datasets with gene symbols can be problematic: https://docs.lamin.ai/faq/symbol-mapping"
            )

        self._dataset = data
        self._obs_fields = categoricals or {}
        self._var_field = var_index
        super().__init__(
            dataset=data,
            categoricals=categoricals,
            sources=sources,
            organism=organism,
            exclude=exclude,
            columns_field=var_index,
        )
        self._obs_df_curator = DataFrameCatCurator(
            df=self._adata.obs,
            categoricals=self.categoricals,
            columns=obs_columns,
            verbosity=verbosity,
            organism=None,
            sources=sources,
            exclude=exclude,
        )

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(self, public: bool = False) -> CurateLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, "var_index": self._var_field},
            public=public,
        )

    def _save_from_var_index(
        self,
        validated_only: bool = True,
    ):
        """Save variable records."""
        update_registry(
            values=list(self._adata.var.index),
            field=self.var_index,
            key="var_index",
            validated_only=validated_only,
            organism=self._organism,
            source=self._sources.get("var_index"),
            exclude=self._exclude.get("var_index"),
        )

    def add_new_from(self, key: str, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records
        """
        self._obs_df_curator.add_new_from(key, **kwargs)

    def add_new_from_var_index(self, **kwargs):
        """Update variable records.

        Args:
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self._save_from_var_index(validated_only=False, **kwargs)

    def validate(self) -> bool:
        """Validate categories.

        This method also registers the validated records in the current instance.

        Args:
            organism: The organism name.

        Returns:
            Whether the AnnData object is validated.
        """
        self._validate_category_error_messages = ""  # reset the error messages

        # add all validated records to the current instance
        self._save_from_var_index(validated_only=True)
        validated_var, non_validated_var = validate_categories(
            self._adata.var.index,
            field=self._var_field,
            key="var_index",
            source=self._sources.get("var_index"),
            hint_print=".add_new_from_var_index()",
            exclude=self._exclude.get("var_index"),
            organism=self._organism,  # type: ignore
        )
        validated_obs = self._obs_df_curator.validate()
        self._non_validated = self._obs_df_curator._non_validated  # type: ignore
        if len(non_validated_var) > 0:
            self._non_validated["var_index"] = non_validated_var  # type: ignore
        self._is_validated = validated_var and validated_obs
        return self._is_validated

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
            self._obs_df_curator.standardize(key)
        # in addition to the obs columns, standardize the var.index
        if key == "var_index" or key == "all":
            syn_mapper = standardize_categories(
                self._adata.var.index,
                field=self.var_index,
                source=self._sources.get("var_index"),
                organism=self._organism,
            )
            if "var_index" in self._non_validated:  # type: ignore
                self._adata.var.index = self._replace_synonyms(
                    "var_index", syn_mapper, self._adata.var.index
                )


class MuDataCatCurator(CatCurator):
    """Curation flow for a ``MuData`` object.

    See also :class:`~lamindb.Curator`.

    Note that if genes or other measurements are removed from the MuData object,
    the object should be recreated using :meth:`~lamindb.Curator.from_mudata`.

    Args:
        mdata: The MuData object to curate.
        var_index: The registry field for mapping the ``.var`` index for each modality.
            For example:
            ``{"modality_1": bt.Gene.ensembl_gene_id, "modality_2": CellMarker.name}``
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
            Use modality keys to specify categoricals for MuData slots such as `"rna:cell_type": bt.CellType.name"`.
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
        ...         "adt": CellMarker.name
        ...     },
        ...     categoricals={
        ...         "cell_type_ontology_id": bt.CellType.ontology_id,
        ...         "donor_id": ULabel.name
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        mdata: MuData,
        var_index: dict[str, FieldAttr],
        categoricals: dict[str, FieldAttr] | None = None,
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
        self._dataset = mdata
        self._organism = organism
        self._var_fields = var_index
        self._columns_field = var_index  # this is for consistency with BaseCatCurator
        self._verify_modality(self._var_fields.keys())
        self._obs_fields = self._parse_categoricals(categoricals)
        self._modalities = set(self._var_fields.keys()) | set(self._obs_fields.keys())
        self._verbosity = verbosity
        self._obs_df_curator = None
        self._organism = organism
        if "obs" in self._modalities:
            self._obs_df_curator = DataFrameCatCurator(
                df=mdata.obs,
                columns=Feature.name,
                categoricals=self._obs_fields.get("obs", {}),
                verbosity=verbosity,
                sources=self._sources.get("obs"),
                exclude=self._exclude.get("obs"),
                organism=organism,
            )
        self._mod_adata_curators = {
            modality: AnnDataCatCurator(
                data=mdata[modality],
                var_index=var_index.get(modality),
                categoricals=self._obs_fields.get(modality),
                verbosity=verbosity,
                sources=self._sources.get(modality),
                exclude=self._exclude.get(modality),
                organism=organism,
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
    def non_validated(self) -> dict[str, dict[str, list[str]]]:  # type: ignore
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return self._non_validated

    def _verify_modality(self, modalities: Iterable[str]):
        """Verify the modality exists."""
        for modality in modalities:
            if modality not in self._dataset.mod.keys():
                raise ValidationError(f"modality '{modality}' does not exist!")

    def _parse_categoricals(self, categoricals: dict[str, FieldAttr]) -> dict:
        """Parse the categorical fields."""
        prefixes = {f"{k}:" for k in self._dataset.mod.keys()}
        obs_fields: dict[str, dict[str, FieldAttr]] = {}
        for k, v in categoricals.items():
            if k not in self._dataset.obs.columns:
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

    def lookup(self, public: bool = False) -> CurateLookup:
        """Lookup categories.

        Args:
            public: Perform lookup on public source ontologies.
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
            public=public,
        )

    @deprecated(new_name="is run by default")
    def add_new_from_columns(
        self,
        modality: str,
        column_names: list[str] | None = None,
        **kwargs,
    ):
        pass

    def add_new_from_var_index(self, modality: str, **kwargs):
        """Update variable records.

        Args:
            modality: The modality name.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self._mod_adata_curators[modality].add_new_from_var_index(**kwargs)

    def _update_registry_all(self):
        """Update all registries."""
        if self._obs_df_curator is not None:
            self._obs_df_curator._update_registry_all(validated_only=True)
        for _, adata_curator in self._mod_adata_curators.items():
            adata_curator._obs_df_curator._update_registry_all(validated_only=True)

    def add_new_from(
        self,
        key: str,
        modality: str | None = None,
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
        modality = modality or "obs"
        if modality in self._mod_adata_curators:
            adata_curator = self._mod_adata_curators[modality]
            adata_curator.add_new_from(key=key, **kwargs)
        if modality == "obs":
            self._obs_df_curator.add_new_from(key=key, **kwargs)

    def validate(self) -> bool:
        """Validate categories."""
        from lamindb.core._settings import settings

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
            obs_validated &= self._obs_df_curator.validate()
            self._non_validated["obs"] = self._obs_df_curator.non_validated  # type: ignore
            logger.print("")

        mods_validated = True
        for modality, adata_curator in self._mod_adata_curators.items():
            logger.info(f'validating categoricals in modality "{modality}"...')
            mods_validated &= adata_curator.validate()
            if len(adata_curator.non_validated) > 0:
                self._non_validated[modality] = adata_curator.non_validated  # type: ignore
            logger.print("")

        self._is_validated = obs_validated & mods_validated
        return self._is_validated

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


def _maybe_curation_keys_not_present(nonval_keys: list[str], name: str):
    if (n := len(nonval_keys)) > 0:
        s = "s" if n > 1 else ""
        are = "are" if n > 1 else "is"
        raise ValidationError(
            f"key{s} passed to {name} {are} not present: {colors.yellow(_format_values(nonval_keys))}"
        )


class TiledbsomaCatCurator(CatCurator):
    """Curation flow for `tiledbsoma.Experiment`.

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
        ...         "donor_id": ULabel.name
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        experiment_uri: UPathStr | Artifact,
        var_index: dict[str, tuple[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
        exclude: dict[str, str | list[str]] | None = None,
    ):
        self._obs_fields = categoricals or {}
        self._var_fields = var_index
        self._columns_field = obs_columns
        if isinstance(experiment_uri, Artifact):
            self._dataset = experiment_uri.path
            self._artifact = experiment_uri
        else:
            self._dataset = UPath(experiment_uri)
            self._artifact = None
        self._organism = organism
        self._sources = sources or {}
        self._exclude = exclude or {}

        self._is_validated: bool | None = False
        self._non_validated_values: dict[str, list] | None = None
        self._validated_values: dict[str, list] = {}
        # filled by _check_save_keys
        self._n_obs: int | None = None
        self._valid_obs_keys: list[str] | None = None
        self._obs_pa_schema: pa.lib.Schema | None = (
            None  # this is needed to create the obs feature set
        )
        self._valid_var_keys: list[str] | None = None
        self._var_fields_flat: dict[str, FieldAttr] | None = None
        self._check_save_keys()

    # check that the provided keys in var_index and categoricals are available in the store
    # and save features
    def _check_save_keys(self):
        from lamindb.core.storage._tiledbsoma import _open_tiledbsoma

        with _open_tiledbsoma(self._dataset, mode="r") as experiment:
            experiment_obs = experiment.obs
            self._n_obs = len(experiment_obs)
            self._obs_pa_schema = experiment_obs.schema
            valid_obs_keys = [
                k for k in self._obs_pa_schema.names if k != "soma_joinid"
            ]
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
        with _open_tiledbsoma(self._dataset, mode="r") as experiment:
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
                    validated_only=True,
                    organism=organism,
                    source=self._sources.get(var_ms_key),
                    exclude=self._exclude.get(var_ms_key),
                )
                _, non_val = validate_categories(
                    values=var_ms_values,
                    field=field,
                    key=var_ms_key,
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
                    validated_only=True,
                    organism=organism,
                    source=self._sources.get(key),
                    exclude=self._exclude.get(key),
                )
                _, non_val = validate_categories(
                    values=values,
                    field=field,
                    key=key,
                    organism=organism,
                    source=self._sources.get(key),
                    exclude=self._exclude.get(key),
                )
                if len(non_val) > 0:
                    validated = False
                    self._non_validated_values[key] = non_val
                else:
                    self._validated_values[key] = values
        self._is_validated = validated
        return self._is_validated

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

    def add_new_from(self, key: str, **kwargs) -> None:
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
                validated_only=False,
                organism=organism,
                source=self._sources.get(k),
                exclude=self._exclude.get(k),
                **kwargs,
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

    def lookup(self, public: bool = False) -> CurateLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CurateLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, **self._var_fields_flat},
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
                source=self._sources.get(k),
                organism=organism,
            )
            if (n_syn_mapper := len(syn_mapper)) == 0:
                continue

            from lamindb.core.storage._tiledbsoma import _open_tiledbsoma

            with _open_tiledbsoma(self._dataset, mode="r") as experiment:
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
            with _open_tiledbsoma(self._dataset, mode="w") as experiment:
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
        *,
        key: str | None = None,
        description: str | None = None,
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

        if not self._is_validated:
            self.validate()
            if not self._is_validated:
                raise ValidationError("Dataset does not validate. Please curate.")

        if self._artifact is None:
            artifact = Artifact(
                self._dataset,
                description=description,
                key=key,
                revises=revises,
                run=run,
            )
            artifact.n_observations = self._n_obs
            artifact.otype = "tiledbsoma"
            artifact.save()
        else:
            artifact = self._artifact

        _schemas_m2m = {}
        if len(self._obs_fields) > 0:
            organism = check_registry_organism(
                self._columns_field.field.model, self._organism
            ).get("organism")
            empty_dict = {field.name: [] for field in self._obs_pa_schema}  # type: ignore
            mock_df = pa.Table.from_pydict(
                empty_dict, schema=self._obs_pa_schema
            ).to_pandas()
            # in parallel to https://github.com/laminlabs/lamindb/blob/2a1709990b5736b480c6de49c0ada47fafc8b18d/lamindb/core/_feature_manager.py#L549-L554
            _schemas_m2m["obs"] = Schema.from_df(
                df=mock_df,
                field=self._columns_field,
                mute=True,
                organism=organism,
            )
        for ms in self._var_fields:
            var_key, var_field = self._var_fields[ms]
            organism = check_registry_organism(
                var_field.field.model, self._organism
            ).get("organism")
            _schemas_m2m[f"{ms}__var"] = Schema.from_values(
                values=self._validated_values[f"{ms}__{var_key}"],
                field=var_field,
                organism=organism,
                raise_validation_error=False,
            )
        artifact._staged__schemas_m2m = _schemas_m2m

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


class SpatialDataCatCurator:
    """Curation flow for a ``Spatialdata`` object.

    See also :class:`~lamindb.Curator`.

    Note that if genes or other measurements are removed from the SpatialData object,
    the object should be recreated.

    In the following docstring, an accessor refers to either a ``.table`` key or the ``sample_metadata_key``.

    Args:
        sdata: The SpatialData object to curate.
        var_index: A dictionary mapping table keys to the ``.var`` indices.
        categoricals: A nested dictionary mapping an accessor to dictionaries that map columns to a registry field.

        organism: The organism name.
        sources: A dictionary mapping an accessor to dictionaries that map columns to Source records.
        exclude: A dictionary mapping an accessor to dictionaries of column names to values to exclude from validation.
            When specific :class:`~bionty.Source` instances are pinned and may lack default values (e.g., "unknown" or "na"),
            using the exclude parameter ensures they are not validated.
        verbosity: The verbosity level of the logger.
        sample_metadata_key: The key in ``.attrs`` that stores the sample level metadata.

    Examples:
        >>> import bionty as bt
        >>> curator = SpatialDataCatCurator(
        ...     sdata,
        ...     var_index={
        ...         "table_1": bt.Gene.ensembl_gene_id,
        ...     },
        ...     categoricals={
        ...         "table1":
        ...             {"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ULabel.name},
        ...         "sample":
        ...             {"experimental_factor": bt.ExperimentalFactor.name},
        ...     },
        ...     organism="human",
        ... )
    """

    def __init__(
        self,
        sdata,
        var_index: dict[str, FieldAttr],
        categoricals: dict[str, dict[str, FieldAttr]] | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, dict[str, Record]] | None = None,
        exclude: dict[str, dict] | None = None,
        *,
        sample_metadata_key: str | None = "sample",
    ) -> None:
        if sources is None:
            sources = {}
        self._sources = sources
        if exclude is None:
            exclude = {}
        self._exclude = exclude
        self._sdata: SpatialData = sdata
        self._sample_metadata_key = sample_metadata_key
        self._organism = organism
        self._var_fields = var_index
        self._verify_accessor_exists(self._var_fields.keys())
        self._categoricals = categoricals
        self._table_keys = set(self._var_fields.keys()) | set(
            self._categoricals.keys() - {self._sample_metadata_key}
        )
        self._verbosity = verbosity
        self._sample_df_curator = None
        if self._sample_metadata_key is not None:
            self._sample_metadata = self._sdata.get_attrs(
                key=self._sample_metadata_key, return_as="df", flatten=True
            )
        self._is_validated = False

        # Check validity of keys in categoricals
        nonval_keys = []
        for accessor, accessor_categoricals in self._categoricals.items():
            if (
                accessor == self._sample_metadata_key
                and self._sample_metadata is not None
            ):
                for key in accessor_categoricals.keys():
                    if key not in self._sample_metadata.columns:
                        nonval_keys.append(key)
            else:
                for key in accessor_categoricals.keys():
                    if key not in self._sdata[accessor].obs.columns:
                        nonval_keys.append(key)

        _maybe_curation_keys_not_present(nonval_keys, "categoricals")

        # check validity of keys in sources and exclude
        for name, dct in (("sources", self._sources), ("exclude", self._exclude)):
            nonval_keys = []
            for accessor, accessor_sources in dct.items():
                if (
                    accessor == self._sample_metadata_key
                    and self._sample_metadata is not None
                ):
                    columns = self._sample_metadata.columns
                elif accessor != self._sample_metadata_key:
                    columns = self._sdata[accessor].obs.columns
                else:
                    continue
                for key in accessor_sources:
                    if key not in columns:
                        nonval_keys.append(key)
            _maybe_curation_keys_not_present(nonval_keys, name)

        # Set up sample level metadata and table Curator objects
        if (
            self._sample_metadata_key is not None
            and self._sample_metadata_key in self._categoricals
        ):
            self._sample_df_curator = DataFrameCatCurator(
                df=self._sample_metadata,
                columns=Feature.name,
                categoricals=self._categoricals.get(self._sample_metadata_key, {}),
                verbosity=verbosity,
                sources=self._sources.get(self._sample_metadata_key),
                exclude=self._exclude.get(self._sample_metadata_key),
                organism=organism,
            )
        self._table_adata_curators = {
            table: AnnDataCatCurator(
                data=sdata[table],
                var_index=var_index.get(table),
                categoricals=self._categoricals.get(table),
                verbosity=verbosity,
                sources=self._sources.get(table),
                exclude=self._exclude.get(table),
                organism=organism,
            )
            for table in self._table_keys
        }

        self._non_validated = None

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry fields to validate variables indices against."""
        return self._var_fields

    @property
    def categoricals(self) -> dict[str, dict[str, FieldAttr]]:
        """Return the categorical keys and fields to validate against."""
        return self._categoricals

    @property
    def non_validated(self) -> dict[str, dict[str, list[str]]]:
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return self._non_validated

    def _verify_accessor_exists(self, accessors: Iterable[str]) -> None:
        """Verify that the accessors exist (either a valid table or in attrs)."""
        for acc in accessors:
            is_present = False
            try:
                self._sdata.get_attrs(key=acc)
                is_present = True
            except KeyError:
                if acc in self._sdata.tables.keys():
                    is_present = True
            if not is_present:
                raise ValidationError(f"Accessor '{acc}' does not exist!")

    def lookup(self, public: bool = False) -> CurateLookup:
        """Look up categories.

        Args:
            public: Whether the lookup is performed on the public reference.
        """
        cat_values_dict = list(self.categoricals.values())[0]
        return CurateLookup(
            categoricals=cat_values_dict,
            slots={"accessors": cat_values_dict.keys()},
            public=public,
        )

    def _update_registry_all(self) -> None:
        """Saves labels of all features for sample and table metadata."""
        if self._sample_df_curator is not None:
            self._sample_df_curator._update_registry_all(
                validated_only=True,
            )
        for _, adata_curator in self._table_adata_curators.items():
            adata_curator._obs_df_curator._update_registry_all(
                validated_only=True,
            )

    def add_new_from_var_index(self, table: str, **kwargs) -> None:
        """Save new values from ``.var.index`` of table.

        Args:
            table: The table key.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        if self._non_validated is None:
            raise ValidationError("Run .validate() first.")
        self._table_adata_curators[table].add_new_from_var_index(**kwargs)
        if table in self.non_validated.keys():
            self._non_validated[table].pop("var_index")

            if len(self.non_validated[table].values()) == 0:
                self.non_validated.pop(table)

    def add_new_from(
        self,
        key: str,
        accessor: str | None = None,
        **kwargs,
    ) -> None:
        """Save new values of categorical from sample level metadata or table.

        Args:
            key: The key referencing the slot in the DataFrame.
            accessor: The accessor key such as 'sample' or 'table x'.
            organism: The organism name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        if self._non_validated is None:
            raise ValidationError("Run .validate() first.")

        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")

        if accessor not in self.categoricals:
            raise ValueError(
                f"Accessor {accessor} is not in 'categoricals'. Include it when creating the SpatialDataCatCurator."
            )

        if accessor in self._table_adata_curators:
            adata_curator = self._table_adata_curators[accessor]
            adata_curator.add_new_from(key=key, **kwargs)
        if accessor == self._sample_metadata_key:
            self._sample_df_curator.add_new_from(key=key, **kwargs)

        if accessor in self.non_validated.keys():
            if len(self.non_validated[accessor].values()) == 0:
                self.non_validated.pop(accessor)

    def standardize(self, key: str, accessor: str | None = None) -> None:
        """Replace synonyms with canonical values.

        Modifies the dataset inplace.

        Args:
            key: The key referencing the slot in the table or sample metadata.
            accessor: The accessor key such as 'sample_key' or 'table_key'.
        """
        if len(self.non_validated) == 0:
            logger.warning("values are already standardized")
            return

        if accessor == self._sample_metadata_key:
            if key not in self._sample_metadata.columns:
                raise ValueError(f"key '{key}' not present in '{accessor}'!")
        else:
            if key not in self._sdata.tables[accessor].obs.columns:
                raise ValueError(f"key '{key}' not present in '{accessor}'!")

        if accessor in self._table_adata_curators.keys():
            adata_curator = self._table_adata_curators[accessor]
            adata_curator.standardize(key)
        if accessor == self._sample_metadata_key:
            self._sample_df_curator.standardize(key)

        if len(self.non_validated[accessor].values()) == 0:
            self.non_validated.pop(accessor)

    def validate(self) -> bool:
        """Validate variables and categorical observations.

        This method also registers the validated records in the current instance:
        - from public sources

        Args:
            organism: The organism name.

        Returns:
            Whether the SpatialData object is validated.
        """
        from lamindb.core._settings import settings

        # add all validated records to the current instance
        verbosity = settings.verbosity
        try:
            settings.verbosity = "error"
            self._update_registry_all()
        finally:
            settings.verbosity = verbosity

        self._non_validated = {}  # type: ignore

        sample_validated = True
        if self._sample_df_curator:
            logger.info(f"validating categoricals of '{self._sample_metadata_key}' ...")
            sample_validated &= self._sample_df_curator.validate()
            if len(self._sample_df_curator.non_validated) > 0:
                self._non_validated["sample"] = self._sample_df_curator.non_validated  # type: ignore
            logger.print("")

        mods_validated = True
        for table, adata_curator in self._table_adata_curators.items():
            logger.info(f"validating categoricals of table '{table}' ...")
            mods_validated &= adata_curator.validate()
            if len(adata_curator.non_validated) > 0:
                self._non_validated[table] = adata_curator.non_validated  # type: ignore
            logger.print("")

        self._is_validated = sample_validated & mods_validated
        return self._is_validated

    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        if not self._is_validated:
            self.validate()
            if not self._is_validated:
                raise ValidationError("Dataset does not validate. Please curate.")

        verbosity = settings.verbosity
        try:
            settings.verbosity = "warning"

            # Write the SpatialData object to a random path in tmp directory
            # The Artifact constructor will move it to the cache
            write_path = f"{settings.cache_dir}/{random.randint(10**7, 10**8 - 1)}.zarr"
            self._sdata.write(write_path)

            # Create the Artifact and associate Artifact metadata
            self._artifact = Artifact(
                write_path,
                description=description,
                key=key,
                revises=revises,
                run=run,
            )
            # According to Tim it is not easy to calculate the number of observations.
            # We would have to write custom code to iterate over labels (which might not even exist at that point)
            self._artifact.otype = "spatialdata"
            self._artifact.save()

            # Link schemas
            feature_kwargs = check_registry_organism(
                (list(self._var_fields.values())[0].field.model),
                self._organism,
            )

            def _add_set_from_spatialdata(
                host: Artifact | Collection | Run,
                var_fields: dict[str, FieldAttr],
                obs_fields: dict[str, FieldAttr] = None,
                mute: bool = False,
                organism: str | Record | None = None,
            ):
                """Add Schemas from SpatialData."""
                if obs_fields is None:
                    obs_fields = {}
                assert host.otype == "spatialdata"  # noqa: S101

                _schemas_m2m = {}

                # sample features
                sample_features = Feature.from_values(self._sample_metadata.columns)  # type: ignore
                if len(sample_features) > 0:
                    _schemas_m2m[self._sample_metadata_key] = Schema(
                        features=sample_features
                    )

                # table features
                for table, field in var_fields.items():
                    table_fs = parse_staged__schemas_m2m_from_anndata(
                        self._sdata[table],
                        var_field=field,
                        obs_field=obs_fields.get(table, Feature.name),
                        mute=mute,
                        organism=organism,
                    )
                    for k, v in table_fs.items():
                        _schemas_m2m[f"['{table}'].{k}"] = v

                def _unify_staged__schemas_m2m_by_hash(
                    _schemas_m2m: MutableMapping[str, Schema],
                ):
                    unique_values: dict[str, Any] = {}

                    for key, value in _schemas_m2m.items():
                        value_hash = (
                            value.hash
                        )  # Assuming each value has a .hash attribute
                        if value_hash in unique_values:
                            _schemas_m2m[key] = unique_values[value_hash]
                        else:
                            unique_values[value_hash] = value

                    return _schemas_m2m

                # link feature sets
                host._staged__schemas_m2m = _unify_staged__schemas_m2m_by_hash(
                    _schemas_m2m
                )
                host.save()

            _add_set_from_spatialdata(
                self._artifact, var_fields=self._var_fields, **feature_kwargs
            )

            # Link labels
            def _add_labels_from_spatialdata(
                data,
                artifact: Artifact,
                fields: dict[str, FieldAttr],
                feature_ref_is_name: bool | None = None,
            ):
                """Add Labels from SpatialData."""
                features = Feature.lookup().dict()
                for key, field in fields.items():
                    feature = features.get(key)
                    registry = field.field.model
                    filter_kwargs = check_registry_organism(registry, self._organism)
                    filter_kwargs_current = get_current_filter_kwargs(
                        registry, filter_kwargs
                    )
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

            for accessor, accessor_fields in self._categoricals.items():
                column_field = self._var_fields.get(accessor)
                if accessor == self._sample_metadata_key:
                    _add_labels_from_spatialdata(
                        self._sample_metadata,
                        self._artifact,
                        accessor_fields,
                        feature_ref_is_name=(
                            None if column_field is None else _ref_is_name(column_field)
                        ),
                    )
                else:
                    _add_labels_from_spatialdata(
                        self._sdata.tables[accessor],
                        self._artifact,
                        accessor_fields,
                        feature_ref_is_name=(
                            None if column_field is None else _ref_is_name(column_field)
                        ),
                    )

        finally:
            settings.verbosity = verbosity

        slug = ln_setup.settings.instance.slug
        if ln_setup.settings.instance.is_remote:  # pragma: no cover
            logger.important(
                f"go to https://lamin.ai/{slug}/artifact/{self._artifact.uid}"
            )

        return self._artifact


def _restrict_obs_fields(
    obs: pd.DataFrame, obs_fields: dict[str, FieldAttr]
) -> dict[str, str]:
    """Restrict the obs fields to name return only available obs fields.

    To simplify the curation, we only validate against either name or ontology_id.
    If both are available, we validate against ontology_id.
    If none are available, we validate against name.
    """
    obs_fields_unique = {k: v for k, v in obs_fields.items() if k in obs.columns}
    for name, field in obs_fields.items():
        if name.endswith("_ontology_term_id"):
            continue
        # if both the ontology id and the name are present, only validate on the ontology_id
        if name in obs.columns and f"{name}_ontology_term_id" in obs.columns:
            obs_fields_unique.pop(name)
        # if the neither name nor ontology id are present, validate on the name
        # this will raise error downstream, we just use name to be more readable
        if name not in obs.columns and f"{name}_ontology_term_id" not in obs.columns:
            obs_fields_unique[name] = field

    # Only retain obs_fields_unique that have keys in adata.obs.columns
    available_obs_fields = {
        k: v for k, v in obs_fields_unique.items() if k in obs.columns
    }

    return available_obs_fields


def _add_defaults_to_obs(
    obs: pd.DataFrame,
    defaults: dict[str, str],
) -> None:
    """Add default columns and values to obs DataFrame."""
    added_defaults: dict = {}
    for name, default in defaults.items():
        if name not in obs.columns and f"{name}_ontology_term_id" not in obs.columns:
            obs[name] = default
            added_defaults[name] = default
            logger.important(
                f"added default value '{default}' to the adata.obs['{name}']"
            )


class CellxGeneAnnDataCatCurator(AnnDataCatCurator):
    """Annotation flow of AnnData based on CELLxGENE schema."""

    _controls_were_created: bool | None = None

    def __init__(
        self,
        adata: ad.AnnData | UPathStr,
        categoricals: dict[str, FieldAttr] | None = None,
        organism: Literal["human", "mouse"] = "human",
        *,
        defaults: dict[str, str] = None,
        extra_sources: dict[str, Record] = None,
        schema_version: Literal["4.0.0", "5.0.0", "5.1.0"] = "5.1.0",
        verbosity: str = "hint",
    ) -> None:
        """CELLxGENE schema curator.

        Args:
            adata: Path to or AnnData object to curate against the CELLxGENE schema.
            categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
                The CELLxGENE Curator maps against the required CELLxGENE fields by default.
            organism: The organism name. CELLxGENE restricts it to 'human' and 'mouse'.
            defaults: Default values that are set if columns or column values are missing.
            extra_sources: A dictionary mapping ``.obs.columns`` to Source records.
                These extra sources are joined with the CELLxGENE fixed sources.
                Use this parameter when subclassing.
            exclude: A dictionary mapping column names to values to exclude.
            schema_version: The CELLxGENE schema version to curate against.
            verbosity: The verbosity level.

        """
        import bionty as bt

        CellxGeneAnnDataCatCurator._init_categoricals_additional_values()

        var_index: FieldAttr = bt.Gene.ensembl_gene_id

        if categoricals is None:
            categoricals = CellxGeneAnnDataCatCurator._get_categoricals()

        self.organism = organism

        VALID_SCHEMA_VERSIONS = {"4.0.0", "5.0.0", "5.1.0"}
        if schema_version not in VALID_SCHEMA_VERSIONS:
            valid_versions = ", ".join(sorted(VALID_SCHEMA_VERSIONS))
            raise ValueError(
                f"Invalid schema_version: {schema_version}. "
                f"Valid versions are: {valid_versions}"
            )
        self.schema_version = schema_version
        self.schema_reference = f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{schema_version}/schema.md"
        with resources.path(
            "lamindb.curators._cellxgene_schemas", "schema_versions.yml"
        ) as schema_versions_path:
            self._pinned_ontologies = _read_schema_versions(schema_versions_path)[
                self.schema_version
            ]

        # Fetch AnnData obs to be able to set defaults and get sources
        if isinstance(adata, ad.AnnData):
            self._adata_obs = adata.obs
        else:
            self._adata_obs = backed_access(upath.create_path(adata)).obs  # type: ignore

        # Add defaults first to ensure that we fetch valid sources
        if defaults:
            _add_defaults_to_obs(self._adata_obs, defaults)

        self.sources = self._create_sources(self._adata_obs)
        self.sources = {
            entity: source
            for entity, source in self.sources.items()
            if source is not None
        }

        # These sources are not a part of the cellxgene schema but rather passed through.
        # This is useful when other Curators extend the CELLxGENE curator
        if extra_sources:
            self.sources = self.sources | extra_sources

        # Exclude default values from validation because they are not available in the pinned sources
        exclude_keys = {
            entity: default
            for entity, default in CellxGeneAnnDataCatCurator._get_categoricals_defaults().items()
            if entity in self._adata_obs.columns  # type: ignore
        }

        super().__init__(
            data=adata,
            var_index=var_index,
            categoricals=_restrict_obs_fields(self._adata_obs, categoricals),
            verbosity=verbosity,
            organism=organism,
            sources=self.sources,
            exclude=exclude_keys,
        )

    @classmethod
    def _init_categoricals_additional_values(cls) -> None:
        import bionty as bt

        import lamindb as ln

        # Note: if you add another control below, be mindful to change the if condition that
        # triggers whether creating these records is re-considered
        if cls._controls_were_created is None:
            cls._controls_were_created = (
                ln.ULabel.filter(name="SuspensionType", is_type=True).one_or_none()
                is not None
            )
        if not cls._controls_were_created:
            logger.important("Creating control labels in the CellxGene schema.")
            bt.CellType(
                ontology_id="unknown",
                name="unknown",
                description="From CellxGene schema.",
            ).save()
            pato = bt.Source.filter(name="pato", version="2024-03-28").one()
            normal = bt.Phenotype.from_source(ontology_id="PATO:0000461", source=pato)
            bt.Disease(
                uid=normal.uid,
                name=normal.name,
                ontology_id=normal.ontology_id,
                description=normal.description,
                source=normal.source,
            ).save()
            bt.Ethnicity(
                ontology_id="na", name="na", description="From CellxGene schema."
            ).save()
            bt.Ethnicity(
                ontology_id="unknown",
                name="unknown",
                description="From CellxGene schema.",
            ).save()
            bt.DevelopmentalStage(
                ontology_id="unknown",
                name="unknown",
                description="From CellxGene schema.",
            ).save()
            bt.Phenotype(
                ontology_id="unknown",
                name="unknown",
                description="From CellxGene schema.",
            ).save()

            tissue_type = ln.ULabel(
                name="TissueType",
                is_type=True,
                description='From CellxGene schema. Is "tissue", "organoid", or "cell culture".',
            ).save()
            ln.ULabel(
                name="tissue", type=tissue_type, description="From CellxGene schema."
            ).save()
            ln.ULabel(
                name="organoid", type=tissue_type, description="From CellxGene schema."
            ).save()
            ln.ULabel(
                name="cell culture",
                type=tissue_type,
                description="From CellxGene schema.",
            ).save()

            suspension_type = ln.ULabel(
                name="SuspensionType",
                is_type=True,
                description='From CellxGene schema. This MUST be "cell", "nucleus", or "na".',
            ).save()
            ln.ULabel(
                name="cell", type=suspension_type, description="From CellxGene schema."
            ).save()
            ln.ULabel(
                name="nucleus",
                type=suspension_type,
                description="From CellxGene schema.",
            ).save()
            ln.ULabel(name="na", type=suspension_type).save()

    @classmethod
    def _get_categoricals(cls) -> dict[str, FieldAttr]:
        import bionty as bt

        return {
            "assay": bt.ExperimentalFactor.name,
            "assay_ontology_term_id": bt.ExperimentalFactor.ontology_id,
            "cell_type": bt.CellType.name,
            "cell_type_ontology_term_id": bt.CellType.ontology_id,
            "development_stage": bt.DevelopmentalStage.name,
            "development_stage_ontology_term_id": bt.DevelopmentalStage.ontology_id,
            "disease": bt.Disease.name,
            "disease_ontology_term_id": bt.Disease.ontology_id,
            # "donor_id": "str",  via pandera
            "self_reported_ethnicity": bt.Ethnicity.name,
            "self_reported_ethnicity_ontology_term_id": bt.Ethnicity.ontology_id,
            "sex": bt.Phenotype.name,
            "sex_ontology_term_id": bt.Phenotype.ontology_id,
            "suspension_type": ULabel.name,
            "tissue": bt.Tissue.name,
            "tissue_ontology_term_id": bt.Tissue.ontology_id,
            "tissue_type": ULabel.name,
            "organism": bt.Organism.name,
            "organism_ontology_term_id": bt.Organism.ontology_id,
        }

    @classmethod
    def _get_categoricals_defaults(cls) -> dict[str, str]:
        return {
            "cell_type": "unknown",
            "development_stage": "unknown",
            "disease": "normal",
            "donor_id": "unknown",
            "self_reported_ethnicity": "unknown",
            "sex": "unknown",
            "suspension_type": "cell",
            "tissue_type": "tissue",
        }

    @property
    def pinned_ontologies(self) -> pd.DataFrame:
        return self._pinned_ontologies

    @property
    def adata(self) -> AnnData:
        return self._adata

    def _create_sources(self, obs: pd.DataFrame) -> dict[str, Record]:
        """Creates a sources dictionary that can be passed to AnnDataCatCurator."""
        import bionty as bt

        # fmt: off
        def _fetch_bionty_source(
            entity: str, organism: str, source: str
        ) -> bt.Source | None:
            """Fetch the Bionty source of the pinned ontology.

            Returns None if the source does not exist.
            """
            version = self._pinned_ontologies.loc[(self._pinned_ontologies.index == entity) &
                                                  (self._pinned_ontologies["organism"] == organism) &
                                                  (self._pinned_ontologies["source"] == source), "version"].iloc[0]
            return bt.Source.filter(organism=organism, entity=f"bionty.{entity}", version=version).first()

        entity_mapping = {
             "var_index": ("Gene", self.organism, "ensembl"),
             "cell_type": ("CellType", "all", "cl"),
             "assay": ("ExperimentalFactor", "all", "efo"),
             "self_reported_ethnicity": ("Ethnicity", self.organism, "hancestro"),
             "development_stage": ("DevelopmentalStage", self.organism, "hsapdv" if self.organism == "human" else "mmusdv"),
             "disease": ("Disease", "all", "mondo"),
             # "organism": ("Organism", "vertebrates", "ensembl"),
             "sex": ("Phenotype", "all", "pato"),
             "tissue": ("Tissue", "all", "uberon"),
        }
        # fmt: on

        # Retain var_index and one of 'entity'/'entity_ontology_term_id' that is present in obs
        entity_to_sources = {
            entity: _fetch_bionty_source(*params)
            for entity, params in entity_mapping.items()
            if entity in obs.columns
            or (f"{entity}_ontology_term_id" in obs.columns and entity != "var_index")
            or entity == "var_index"
        }

        return entity_to_sources

    def _convert_name_to_ontology_id(self, values: pd.Series, field: FieldAttr):
        """Converts a column that stores a name into a column that stores the ontology id.

        cellxgene expects the obs columns to be {entity}_ontology_id columns and disallows {entity} columns.
        """
        field_name = field.field.name
        assert field_name == "name"  # noqa: S101
        cols = ["name", "ontology_id"]
        registry = field.field.model

        if hasattr(registry, "ontology_id"):
            validated_records = registry.filter(**{f"{field_name}__in": values})
            mapper = (
                pd.DataFrame(validated_records.values_list(*cols))
                .set_index(0)
                .to_dict()[1]
            )
            return values.map(mapper)

    def validate(self) -> bool:  # type: ignore
        """Validates the AnnData object against most cellxgene requirements."""
        # Verify that all required obs columns are present
        missing_obs_fields = [
            name
            for name in CellxGeneAnnDataCatCurator._get_categoricals_defaults().keys()
            if name not in self._adata.obs.columns
            and f"{name}_ontology_term_id" not in self._adata.obs.columns
        ]
        if len(missing_obs_fields) > 0:
            missing_obs_fields_str = ", ".join(list(missing_obs_fields))
            logger.error(f"missing required obs columns {missing_obs_fields_str}")
            logger.info(
                "consider initializing a Curate object like 'Curate(adata, defaults=cxg.CellxGeneAnnDataCatCurator._get_categoricals_defaults())'"
                "to automatically add these columns with default values."
            )
            return False

        # Verify that no cellxgene reserved names are present
        reserved_names = {
            "ethnicity",
            "ethnicity_ontology_term_id",
            "X_normalization",
            "default_field",
            "layer_descriptions",
            "tags",
            "versions",
            "contributors",
            "preprint_doi",
            "project_description",
            "project_links",
            "project_name",
            "publication_doi",
        }
        matched_columns = [
            column for column in self._adata.obs.columns if column in reserved_names
        ]
        if len(matched_columns) > 0:
            raise ValueError(
                f"AnnData object must not contain obs columns {matched_columns} which are"
                " reserved from previous schema versions."
            )

        return super().validate()

    def to_cellxgene_anndata(
        self, is_primary_data: bool, title: str | None = None
    ) -> ad.AnnData:
        """Converts the AnnData object to the cellxgene-schema input format.

        cellxgene expects the obs fields to be {entity}_ontology_id fields and has many further requirements which are
        documented here: https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema.
        This function checks for most but not all requirements of the CELLxGENE schema.
        If you want to ensure that it fully adheres to the CELLxGENE schema, run `cellxgene-schema` on the AnnData object.

        Args:
            is_primary_data: Whether the measured data is primary data or not.
            title: Title of the AnnData object. Commonly the name of the publication.

        Returns:
            An AnnData object which adheres to the cellxgene-schema.
        """
        # Create a copy since we modify the AnnData object extensively
        adata_cxg = self._adata.copy()

        # cellxgene requires an embedding
        embedding_pattern = r"^[a-zA-Z][a-zA-Z0-9_.-]*$"
        exclude_key = "spatial"
        matching_keys = [
            key
            for key in adata_cxg.obsm.keys()
            if re.match(embedding_pattern, key) and key != exclude_key
        ]
        if len(matching_keys) == 0:
            raise ValueError(
                "Unable to find an embedding key. Please calculate an embedding."
            )

        # convert name column to ontology_term_id column
        for column in adata_cxg.obs.columns:
            if column in self.categoricals and not column.endswith("_ontology_term_id"):
                mapped_column = self._convert_name_to_ontology_id(
                    adata_cxg.obs[column], field=self.categoricals.get(column)
                )
                if mapped_column is not None:
                    adata_cxg.obs[f"{column}_ontology_term_id"] = mapped_column

        # drop the name columns for ontologies. cellxgene does not allow them.
        drop_columns = [
            i
            for i in adata_cxg.obs.columns
            if f"{i}_ontology_term_id" in adata_cxg.obs.columns
        ]
        adata_cxg.obs.drop(columns=drop_columns, inplace=True)

        # Add cellxgene metadata to AnnData object
        if "is_primary_data" not in adata_cxg.obs.columns:
            adata_cxg.obs["is_primary_data"] = is_primary_data
        if "feature_is_filtered" not in adata_cxg.var.columns:
            logger.warn(
                "column 'feature_is_filtered' not present in var. Setting to default"
                " value of False."
            )
            adata_cxg.var["feature_is_filtered"] = False
        if title is None:
            raise ValueError("please pass a title!")
        else:
            adata_cxg.uns["title"] = title
        adata_cxg.uns["cxg_lamin_schema_reference"] = self.schema_reference
        adata_cxg.uns["cxg_lamin_schema_version"] = self.schema_version

        return adata_cxg


class ValueUnit:
    """Base class for handling value-unit combinations."""

    @staticmethod
    def parse_value_unit(value: str, is_dose: bool = True) -> tuple[str, str] | None:
        """Parse a string containing a value and unit into a tuple."""
        if not isinstance(value, str) or not value.strip():
            return None

        value = str(value).strip()
        match = re.match(r"^(\d*\.?\d{0,1})\s*([a-zA-ZμµΜ]+)$", value)

        if not match:
            raise ValueError(
                f"Invalid format: {value}. Expected format: number with max 1 decimal place + unit"
            )

        number, unit = match.groups()
        formatted_number = f"{float(number):.1f}"

        if is_dose:
            standardized_unit = DoseHandler.standardize_unit(unit)
            if not DoseHandler.validate_unit(standardized_unit):
                raise ValueError(
                    f"Invalid dose unit: {unit}. Must be convertible to one of: nM, μM, mM, M"
                )
        else:
            standardized_unit = TimeHandler.standardize_unit(unit)
            if not TimeHandler.validate_unit(standardized_unit):
                raise ValueError(
                    f"Invalid time unit: {unit}. Must be convertible to one of: h, m, s, d, y"
                )

        return formatted_number, standardized_unit


class DoseHandler:
    """Handler for dose-related operations."""

    VALID_UNITS = {"nM", "μM", "µM", "mM", "M"}
    UNIT_MAP = {
        "nm": "nM",
        "NM": "nM",
        "um": "μM",
        "UM": "μM",
        "μm": "μM",
        "μM": "μM",
        "µm": "μM",
        "µM": "μM",
        "mm": "mM",
        "MM": "mM",
        "m": "M",
        "M": "M",
    }

    @classmethod
    def validate_unit(cls, unit: str) -> bool:
        """Validate if the dose unit is acceptable."""
        return unit in cls.VALID_UNITS

    @classmethod
    def standardize_unit(cls, unit: str) -> str:
        """Standardize dose unit to standard formats."""
        return cls.UNIT_MAP.get(unit, unit)

    @classmethod
    def validate_values(cls, values: pd.Series) -> list:
        """Validate pert_dose values with strict case checking."""
        errors = []

        for idx, value in values.items():
            if pd.isna(value):
                continue

            if isinstance(value, (int, float)):
                errors.append(
                    f"Row {idx} - Missing unit for dose: {value}. Must include a unit (nM, μM, mM, M)"
                )
                continue

            try:
                ValueUnit.parse_value_unit(value, is_dose=True)
            except ValueError as e:
                errors.append(f"Row {idx} - {str(e)}")

        return errors


class TimeHandler:
    """Handler for time-related operations."""

    VALID_UNITS = {"h", "m", "s", "d", "y"}

    @classmethod
    def validate_unit(cls, unit: str) -> bool:
        """Validate if the time unit is acceptable."""
        return unit == unit.lower() and unit in cls.VALID_UNITS

    @classmethod
    def standardize_unit(cls, unit: str) -> str:
        """Standardize time unit to standard formats."""
        if unit.startswith("hr"):
            return "h"
        elif unit.startswith("min"):
            return "m"
        elif unit.startswith("sec"):
            return "s"
        return unit[0].lower()

    @classmethod
    def validate_values(cls, values: pd.Series) -> list:
        """Validate pert_time values."""
        errors = []

        for idx, value in values.items():
            if pd.isna(value):
                continue

            if isinstance(value, (int, float)):
                errors.append(
                    f"Row {idx} - Missing unit for time: {value}. Must include a unit (h, m, s, d, y)"
                )
                continue

            try:
                ValueUnit.parse_value_unit(value, is_dose=False)
            except ValueError as e:
                errors.append(f"Row {idx} - {str(e)}")

        return errors


class PertAnnDataCatCurator(CellxGeneAnnDataCatCurator):
    """Curator flow for Perturbation data."""

    PERT_COLUMNS = {"compound", "genetic", "biologic", "physical"}

    def __init__(
        self,
        adata: ad.AnnData,
        organism: Literal["human", "mouse"] = "human",
        pert_dose: bool = True,
        pert_time: bool = True,
        *,
        verbosity: str = "hint",
        cxg_schema_version: Literal["5.0.0", "5.1.0"] = "5.1.0",
    ):
        """Initialize the curator with configuration and validation settings."""
        import bionty as bt

        self._pert_time = pert_time
        self._pert_dose = pert_dose

        self._validate_initial_data(adata)
        self._setup_configuration(adata)

        self._setup_sources(adata)
        self._setup_compound_source()

        super().__init__(
            adata=adata,
            categoricals=self.PT_CATEGORICALS,
            defaults=self.PT_DEFAULT_VALUES,
            verbosity=verbosity,
            organism=organism,
            extra_sources=self.PT_SOURCES,
            schema_version=cxg_schema_version,
        )

    def _setup_configuration(self, adata: ad.AnnData):
        """Set up default configuration values."""
        import bionty as bt
        import wetlab as wl

        self.PT_DEFAULT_VALUES = (
            CellxGeneAnnDataCatCurator._get_categoricals_defaults()
            | {
                "cell_line": "unknown",
                "pert_target": "unknown",
            }
        )

        self.PT_CATEGORICALS = CellxGeneAnnDataCatCurator._get_categoricals() | {
            k: v
            for k, v in {
                "cell_line": bt.CellLine.name,
                "pert_target": wl.PerturbationTarget.name,
                "pert_genetic": wl.GeneticPerturbation.name,
                "pert_compound": wl.Compound.name,
                "pert_biologic": wl.Biologic.name,
                "pert_physical": wl.EnvironmentalPerturbation.name,
            }.items()
            if k in adata.obs.columns
        }
        # if "donor_id" in self.PT_CATEGORICALS:
        #     self.PT_CATEGORICALS["donor_id"] = Donor.name

    def _setup_sources(self, adata: ad.AnnData):
        """Set up data sources."""
        self.PT_SOURCES = {}
        # if "cell_line" in adata.obs.columns:
        #     self.PT_SOURCES["cell_line"] = (
        #         bt.Source.filter(name="depmap").first()
        #     )
        if "pert_compound" in adata.obs.columns:
            import bionty as bt

            self.PT_SOURCES["pert_compound"] = bt.Source.filter(
                entity="wetlab.Compound", name="chebi"
            ).first()

    def _validate_initial_data(self, adata: ad.AnnData):
        """Validate the initial data structure."""
        self._validate_required_columns(adata)
        self._validate_perturbation_types(adata)

    def _validate_required_columns(self, adata: ad.AnnData):
        """Validate required columns are present."""
        if "pert_target" not in adata.obs.columns:
            if (
                "pert_name" not in adata.obs.columns
                or "pert_type" not in adata.obs.columns
            ):
                raise ValidationError(
                    "either 'pert_target' or both 'pert_name' and 'pert_type' must be present"
                )
        else:
            if "pert_name" not in adata.obs.columns:
                logger.warning(
                    "no 'pert' column found in adata.obs, will only curate 'pert_target'"
                )
            elif "pert_type" not in adata.obs.columns:
                raise ValidationError("both 'pert' and 'pert_type' must be present")

    def _validate_perturbation_types(self, adata: ad.AnnData):
        """Validate perturbation types."""
        if "pert_type" in adata.obs.columns:
            data_pert_types = set(adata.obs["pert_type"].unique())
            invalid_pert_types = data_pert_types - self.PERT_COLUMNS
            if invalid_pert_types:
                raise ValidationError(
                    f"invalid pert_type found: {invalid_pert_types}!\n"
                    f"    → allowed values: {self.PERT_COLUMNS}"
                )
            self._process_perturbation_types(adata, data_pert_types)

    def _process_perturbation_types(self, adata: ad.AnnData, pert_types: set):
        """Process and map perturbation types."""
        for pert_type in pert_types:
            col_name = "pert_" + pert_type
            adata.obs[col_name] = adata.obs["pert_name"].where(
                adata.obs["pert_type"] == pert_type, None
            )
            if adata.obs[col_name].dtype.name == "category":
                adata.obs[col_name].cat.remove_unused_categories()
            logger.important(f"mapped 'pert_name' to '{col_name}'")

    def _setup_compound_source(self):
        """Set up the compound source with muted logging."""
        import bionty as bt
        import wetlab as wl

        with logger.mute():
            chebi_source = bt.Source.filter(
                entity="wetlab.Compound", name="chebi"
            ).first()
            if not chebi_source:
                wl.Compound.add_source(
                    bt.Source.filter(entity="Drug", name="chebi").first()
                )

    def validate(self) -> bool:  # type: ignore
        """Validate the AnnData object."""
        validated = super().validate()

        if self._pert_dose:
            validated &= self._validate_dose_column()
        if self._pert_time:
            validated &= self._validate_time_column()

        self._is_validated = validated

        # sort columns
        first_columns = [
            "pert_target",
            "pert_genetic",
            "pert_compound",
            "pert_biologic",
            "pert_physical",
            "pert_dose",
            "pert_time",
            "organism",
            "cell_line",
            "cell_type",
            "disease",
            "tissue_type",
            "tissue",
            "assay",
            "suspension_type",
            "donor_id",
            "sex",
            "self_reported_ethnicity",
            "development_stage",
            "pert_name",
            "pert_type",
        ]
        sorted_columns = [
            col for col in first_columns if col in self._adata.obs.columns
        ] + [col for col in self._adata.obs.columns if col not in first_columns]
        # must assign to self._df to ensure .standardize works correctly
        self._obs_df = self._adata.obs[sorted_columns]
        self._adata.obs = self._obs_df
        return validated

    def standardize(self, key: str) -> pd.DataFrame:
        """Standardize the AnnData object."""
        super().standardize(key)
        self._adata.obs = self._obs_df

    def _validate_dose_column(self) -> bool:
        """Validate the dose column."""
        if not Feature.filter(name="pert_dose").exists():
            Feature(name="pert_dose", dtype="str").save()  # type: ignore

        dose_errors = DoseHandler.validate_values(self._adata.obs["pert_dose"])
        if dose_errors:
            self._log_validation_errors("pert_dose", dose_errors)
            return False
        return True

    def _validate_time_column(self) -> bool:
        """Validate the time column."""
        if not Feature.filter(name="pert_time").exists():
            Feature(name="pert_time", dtype="str").save()  # type: ignore

        time_errors = TimeHandler.validate_values(self._adata.obs["pert_time"])
        if time_errors:
            self._log_validation_errors("pert_time", time_errors)
            return False
        return True

    def _log_validation_errors(self, column: str, errors: list):
        """Log validation errors with formatting."""
        errors_print = "\n    ".join(errors)
        logger.warning(
            f"invalid {column} values found!\n    {errors_print}\n"
            f"    → run {colors.cyan('standardize_dose_time()')}"
        )

    def standardize_dose_time(self) -> pd.DataFrame:
        """Standardize dose and time values."""
        standardized_df = self._adata.obs.copy()

        if "pert_dose" in self._adata.obs.columns:
            standardized_df = self._standardize_column(
                standardized_df, "pert_dose", is_dose=True
            )

        if "pert_time" in self._adata.obs.columns:
            standardized_df = self._standardize_column(
                standardized_df, "pert_time", is_dose=False
            )

        self._adata.obs = standardized_df
        return standardized_df

    def _standardize_column(
        self, df: pd.DataFrame, column: str, is_dose: bool
    ) -> pd.DataFrame:
        """Standardize values in a specific column."""
        for idx, value in self._adata.obs[column].items():
            if pd.isna(value) or (
                isinstance(value, str) and (not value.strip() or value.lower() == "nan")
            ):
                df.at[idx, column] = None
                continue

            try:
                num, unit = ValueUnit.parse_value_unit(value, is_dose=is_dose)
                df.at[idx, column] = f"{num}{unit}"
            except ValueError:
                continue

        return df


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
    organism: str | None = None,
    source: Record | None = None,
    exclude: str | list | None = None,
    hint_print: str | None = None,
    curator: CatCurator | None = None,
) -> tuple[bool, list]:
    """Validate ontology terms in a pandas series using LaminDB registries.

    Args:
        values: The values to validate.
        field: The field attribute.
        key: The key referencing the slot in the DataFrame.
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

    # inspect the non-validated values from public (bionty only)
    values_validated = []
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
        if curator is not None:
            curator._validate_category_error_messages = strip_ansi_codes(
                warning_message
            )
        logger.indent = ""
        return False, non_validated


def standardize_categories(
    values: Iterable[str],
    field: FieldAttr,
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
    return syn_mapper


def validate_categories_in_df(
    df: pd.DataFrame,
    fields: dict[str, FieldAttr],
    sources: dict[str, Record] = None,
    exclude: dict | None = None,
    curator: CatCurator | None = None,
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
            source=sources.get(key),
            exclude=exclude.get(key) if exclude else None,
            curator=curator,
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
    key: str | None = None,
    revises: Artifact | None = None,
    run: Run | None = None,
    schema: Schema | None = None,
) -> Artifact:
    """Save all metadata with an Artifact.

    Args:
        data: The DataFrame/AnnData/MuData object to save.
        fields: A dictionary mapping obs_column to registry_field.
        columns_field: The registry field to validate variables index against.
        description: A description of the artifact.
        organism: The organism name.
        type: The artifact type.
        key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
        revises: Previous version of the artifact. Triggers a revision.
        run: The run that creates the artifact.

    Returns:
        The saved Artifact.
    """
    from .._artifact import data_is_anndata, data_is_mudata
    from ..core._data import add_labels

    print(data)
    if data_is_anndata(data):
        artifact = Artifact.from_anndata(
            data, description=description, key=key, revises=revises, run=run
        )
    elif isinstance(data, pd.DataFrame):
        artifact = Artifact.from_df(
            data, description=description, key=key, revises=revises, run=run
        )
    elif data_is_mudata(data):
        artifact = Artifact.from_mudata(
            data,
            description=description,
            key=key,
            revises=revises,
            run=run,
        )
    artifact.save()

    if organism is not None:
        feature_kwargs = check_registry_organism(
            (
                list(columns_field.values())[0].field.model
                if isinstance(columns_field, dict)
                else columns_field.field.model
            ),
            organism,
        )
    else:
        feature_kwargs = {}

    if artifact.otype == "DataFrame":
        # old style
        artifact.features._add_set_from_df(field=columns_field, **feature_kwargs)  # type: ignore
        # new style (doing both for the time being)
        artifact.schema = schema
    elif artifact.otype == "AnnData":
        artifact.features._add_set_from_anndata(  # type: ignore
            var_field=columns_field, **feature_kwargs
        )
    elif artifact.otype == "MuData":
        artifact.features._add_set_from_mudata(  # type: ignore
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
            # multi-value columns are separated by "|"
            if df[key].str.contains("|").any():
                values = df[key].str.split("|").explode().unique()
            else:
                values = df[key].unique()
            labels = registry.from_values(
                values,
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

    if artifact.otype == "MuData":
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
    if ln_setup.settings.instance.is_remote:  # pdagma: no cover
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
    validated_only: bool = True,
    df: pd.DataFrame | None = None,
    organism: str | None = None,
    dtype: str | None = None,
    source: Record | None = None,
    exclude: str | list | None = None,
    **kwargs,
) -> None:
    """Save features or labels records in the default instance..

    Args:
        values: A list of values to be saved as labels.
        field: The FieldAttr object representing the field for which labels are being saved.
        key: The name of the feature to save.
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
    values = [i for i in values if isinstance(i, str) and i]
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

        # save non-validated/new records
        labels_saved["new"] = non_validated_labels
        if not validated_only:
            non_validated_records: RecordList[Any] = []  # type: ignore
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
    from .._from_values import _format_values

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
    from .._can_curate import get_name_field

    name_field = get_name_field(field.field.model)
    return field.field.name == name_field


# backward compat constructors ------------------


@classmethod  # type: ignore
def from_df(
    cls,
    df: pd.DataFrame,
    categoricals: dict[str, FieldAttr] | None = None,
    columns: FieldAttr = Feature.name,
    verbosity: str = "hint",
    organism: str | None = None,
) -> DataFrameCatCurator:
    return DataFrameCatCurator(
        df=df,
        categoricals=categoricals,
        columns=columns,
        verbosity=verbosity,
        organism=organism,
    )


@classmethod  # type: ignore
def from_anndata(
    cls,
    data: ad.AnnData | UPathStr,
    var_index: FieldAttr,
    categoricals: dict[str, FieldAttr] | None = None,
    obs_columns: FieldAttr = Feature.name,
    verbosity: str = "hint",
    organism: str | None = None,
    sources: dict[str, Record] | None = None,
) -> AnnDataCatCurator:
    return AnnDataCatCurator(
        data=data,
        var_index=var_index,
        categoricals=categoricals,
        obs_columns=obs_columns,
        verbosity=verbosity,
        organism=organism,
        sources=sources,
    )


@classmethod  # type: ignore
def from_mudata(
    cls,
    mdata: MuData,
    var_index: dict[str, dict[str, FieldAttr]],
    categoricals: dict[str, FieldAttr] | None = None,
    verbosity: str = "hint",
    organism: str | None = None,
) -> MuDataCatCurator:
    return MuDataCatCurator(
        mdata=mdata,
        var_index=var_index,
        categoricals=categoricals,
        verbosity=verbosity,
        organism=organism,
    )


@classmethod  # type: ignore
def from_tiledbsoma(
    cls,
    experiment_uri: UPathStr,
    var_index: dict[str, tuple[str, FieldAttr]],
    categoricals: dict[str, FieldAttr] | None = None,
    obs_columns: FieldAttr = Feature.name,
    organism: str | None = None,
    sources: dict[str, Record] | None = None,
    exclude: dict[str, str | list[str]] | None = None,
) -> TiledbsomaCatCurator:
    return TiledbsomaCatCurator(
        experiment_uri=experiment_uri,
        var_index=var_index,
        categoricals=categoricals,
        obs_columns=obs_columns,
        organism=organism,
        sources=sources,
        exclude=exclude,
    )


@classmethod  # type: ignore
def from_spatialdata(
    cls,
    sdata,
    var_index: dict[str, FieldAttr],
    categoricals: dict[str, dict[str, FieldAttr]] | None = None,
    organism: str | None = None,
    sources: dict[str, dict[str, Record]] | None = None,
    exclude: dict[str, dict] | None = None,
    verbosity: str = "hint",
    *,
    sample_metadata_key: str = "sample",
):
    try:
        import spatialdata
    except ImportError as e:
        raise ImportError("Please install spatialdata: pip install spatialdata") from e

    return SpatialDataCatCurator(
        sdata=sdata,
        var_index=var_index,
        categoricals=categoricals,
        verbosity=verbosity,
        organism=organism,
        sources=sources,
        exclude=exclude,
        sample_metadata_key=sample_metadata_key,
    )


Curator.from_df = from_df  # type: ignore
Curator.from_anndata = from_anndata  # type: ignore
Curator.from_mudata = from_mudata  # type: ignore
Curator.from_spatialdata = from_spatialdata  # type: ignore
Curator.from_tiledbsoma = from_tiledbsoma  # type: ignore
