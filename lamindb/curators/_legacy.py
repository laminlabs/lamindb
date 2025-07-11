from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING, Any, Iterable

import pandas as pd
import pyarrow as pa
from lamin_utils import colors, logger
from lamindb_setup.core import deprecated
from lamindb_setup.core.upath import UPath

from lamindb.core._compat import is_package_installed
from lamindb.models.artifact import data_is_scversedatastructure

from ..errors import InvalidArgument

if TYPE_CHECKING:
    from lamindb_setup.types import UPathStr
    from mudata import MuData
    from spatialdata import SpatialData

    from lamindb.models import SQLRecord
from lamindb.base.types import FieldAttr  # noqa
from lamindb.models import (
    Artifact,
    Feature,
    SQLRecord,
    Run,
    Schema,
)
from lamindb.models.artifact import (
    add_labels,
)
from lamindb.models._from_values import _format_values
from .core import CatLookup, CatVector
from ..errors import ValidationError
import anndata as ad


def _ref_is_name(field: FieldAttr | None) -> bool | None:
    """Check if the reference field is a name field."""
    from ..models.can_curate import get_name_field

    if field is not None:
        name_field = get_name_field(field.field.model)
        return field.field.name == name_field
    return None


class CatManager:
    """Manage categoricals by updating registries.

    This class is accessible from within a `DataFrameCurator` via the `.cat` attribute.

    If you find non-validated values, you have several options:

    - new values found in the data can be registered via `DataFrameCurator.cat.add_new_from()` :meth:`~lamindb.curators.DataFrameCatManager.add_new_from`
    - non-validated values can be accessed via `DataFrameCurator.cat.add_new_from()` :meth:`~lamindb.curators.DataFrameCatManager.non_validated` and addressed manually
    """

    def __init__(self, *, dataset, categoricals, sources, columns_field=None):
        # the below is shared with Curator
        self._artifact: Artifact = None  # pass the dataset as an artifact
        self._dataset: Any = dataset  # pass the dataset as a UPathStr or data object
        if isinstance(self._dataset, Artifact):
            self._artifact = self._dataset
            if self._artifact.otype in {"DataFrame", "AnnData"}:
                self._dataset = self._dataset.load(
                    is_run_input=False  # we already track this in the Curator constructor
                )
        self._is_validated: bool = False
        # shared until here
        self._categoricals = categoricals or {}
        self._non_validated = None
        self._sources = sources or {}
        self._columns_field = columns_field
        self._validate_category_error_messages: str = ""
        self._cat_vectors: dict[str, CatVector] = {}

    @property
    def non_validated(self) -> dict[str, list[str]]:
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        return {
            key: cat_vector._non_validated
            for key, cat_vector in self._cat_vectors.items()
            if cat_vector._non_validated and key != "columns"
        }

    @property
    def categoricals(self) -> dict:
        """Return the columns fields to validate against."""
        return self._categoricals

    def validate(self) -> bool:
        """Validate dataset.

        This method also registers the validated records in the current instance.

        Returns:
            The boolean `True` if the dataset is validated. Otherwise, a string with the error message.
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
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
        # Make sure all labels are saved in the current instance
        if not self._is_validated:
            self.validate()  # returns True or False
            if not self._is_validated:  # need to raise error manually
                raise ValidationError("Dataset does not validate. Please curate.")

        if self._artifact is None:
            if isinstance(self._dataset, pd.DataFrame):
                artifact = Artifact.from_df(
                    self._dataset,
                    key=key,
                    description=description,
                    revises=revises,
                    run=run,
                )
            elif isinstance(self._dataset, ad.AnnData):
                artifact = Artifact.from_anndata(
                    self._dataset,
                    key=key,
                    description=description,
                    revises=revises,
                    run=run,
                )
            elif data_is_scversedatastructure(self._dataset, "MuData"):
                artifact = Artifact.from_mudata(
                    self._dataset,
                    key=key,
                    description=description,
                    revises=revises,
                    run=run,
                )
            elif data_is_scversedatastructure(self._dataset, "SpatialData"):
                artifact = Artifact.from_spatialdata(
                    self._dataset,
                    key=key,
                    description=description,
                    revises=revises,
                    run=run,
                )
            else:
                raise InvalidArgument(  # pragma: no cover
                    "data must be one of pd.Dataframe, AnnData, MuData, SpatialData."
                )
            self._artifact = artifact.save()

        legacy_annotate_artifact(  # type: ignore
            self._artifact,
            index_field=self._columns_field,
            cat_vectors=self._cat_vectors,
        )
        return self._artifact


class DataFrameCatManager(CatManager):
    """Categorical manager for `DataFrame`."""

    def __init__(
        self,
        df: pd.DataFrame | Artifact,
        columns_field: FieldAttr = Feature.name,
        columns_names: Iterable[str] | None = None,
        categoricals: dict[str, FieldAttr] | None = None,
        sources: dict[str, SQLRecord] | None = None,
        index: Feature | None = None,
    ) -> None:
        self._non_validated = None
        self._index = index
        super().__init__(
            dataset=df,
            columns_field=columns_field,
            categoricals=categoricals,
            sources=sources,
        )
        if columns_names is None:
            columns_names = []
        if columns_field == Feature.name:
            values = list(self._categoricals.keys())  # backward compat
            self._cat_vectors["columns"] = CatVector(
                values_getter=values,
                field=self._columns_field,
                key="columns" if isinstance(self._dataset, pd.DataFrame) else "keys",
                source=self._sources.get("columns"),
            )
            if isinstance(self._categoricals, dict):  # backward compat
                self._cat_vectors["columns"].validate()
        else:
            # NOTE: for var_index right now
            self._cat_vectors["columns"] = CatVector(
                values_getter=lambda: self._dataset.columns,  # lambda ensures the inplace update
                values_setter=lambda new_values: setattr(
                    self._dataset, "columns", pd.Index(new_values)
                ),
                field=self._columns_field,
                key="columns",
                source=self._sources.get("columns"),
            )
        self._cat_vectors["columns"].add_new()
        for key, field in self._categoricals.items():
            self._cat_vectors[key] = CatVector(
                values_getter=lambda k=key: self._dataset[
                    k
                ],  # Capture key as default argument
                values_setter=lambda new_values, k=key: self._dataset.__setitem__(
                    k, new_values
                ),
                field=field,
                key=key,
                source=self._sources.get(key),
                feature=Feature.get(name=key),
            )

    def lookup(self, public: bool = False) -> CatLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CatLookup(
            categoricals=self._categoricals,
            slots={"columns": self._columns_field},
            public=public,
            sources=self._sources,
        )

    def validate(self) -> bool:
        """Validate variables and categorical observations."""
        self._validate_category_error_messages = ""  # reset the error messages

        validated = True
        for _, cat_vector in self._cat_vectors.items():
            cat_vector.validate()
            validated &= cat_vector.is_validated
        self._is_validated = validated
        self._non_validated = {}  # so it's no longer None

        if self._index is not None:
            # cat_vector.validate() populates validated labels
            # the index should become part of the feature set corresponding to the dataframe
            self._cat_vectors["columns"].labels.insert(0, self._index)  # type: ignore

        return self._is_validated

    def standardize(self, key: str) -> None:
        """Replace synonyms with standardized values.

        Modifies the input dataset inplace.

        Args:
            key: The key referencing the column in the DataFrame to standardize.
        """
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")

        if key == "all":
            logger.warning(
                "'all' is deprecated, please pass a single key from `.non_validated.keys()` instead!"
            )
            for k in self.non_validated.keys():
                self._cat_vectors[k].standardize()
        else:
            self._cat_vectors[key].standardize()

    def add_new_from(self, key: str, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            **kwargs: Additional keyword arguments to pass to create new records
        """
        if len(kwargs) > 0 and key == "all":
            raise ValueError("Cannot pass additional arguments to 'all' key!")
        if key == "all":
            logger.warning(
                "'all' is deprecated, please pass a single key from `.non_validated.keys()` instead!"
            )
            for k in self.non_validated.keys():
                self._cat_vectors[k].add_new(**kwargs)
        else:
            self._cat_vectors[key].add_new(**kwargs)

    @deprecated(
        new_name="Run.filter(transform=context.run.transform, output_artifacts=None)"
    )
    def clean_up_failed_runs(self):
        """Clean up previous failed runs that don't save any outputs."""
        from lamindb.core._context import context

        if context.run is not None:
            Run.filter(transform=context.run.transform, output_artifacts=None).exclude(
                uid=context.run.uid
            ).delete()


class AnnDataCatManager(CatManager):
    """Categorical manager for `AnnData`."""

    def __init__(
        self,
        data: ad.AnnData | Artifact,
        var_index: FieldAttr | None = None,
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        sources: dict[str, SQLRecord] | None = None,
    ) -> None:
        if isinstance(var_index, str):
            raise TypeError(
                "var_index parameter has to be a field, e.g. Gene.ensembl_gene_id"
            )

        if not data_is_scversedatastructure(data, "AnnData"):
            raise TypeError("data has to be an AnnData object")

        if "symbol" in str(var_index):
            logger.warning(
                "indexing datasets with gene symbols can be problematic: https://docs.lamin.ai/faq/symbol-mapping"
            )

        self._obs_fields = categoricals or {}
        self._var_field = var_index
        self._sources = sources or {}
        super().__init__(
            dataset=data,
            categoricals=categoricals,
            sources=self._sources,
            columns_field=var_index,
        )
        self._adata = self._dataset
        self._obs_df_curator = DataFrameCatManager(
            df=self._adata.obs,
            categoricals=self.categoricals,
            columns_field=obs_columns,
            sources=self._sources,
        )
        self._cat_vectors = self._obs_df_curator._cat_vectors.copy()
        if var_index is not None:
            self._cat_vectors["var_index"] = CatVector(
                values_getter=lambda: self._adata.var.index,
                values_setter=lambda new_values: setattr(
                    self._adata.var, "index", pd.Index(new_values)
                ),
                field=self._var_field,
                key="var_index",
                source=self._sources.get("var_index"),
            )

    @property
    def var_index(self) -> FieldAttr:
        """Return the registry field to validate variables index against."""
        return self._var_field

    @property
    def categoricals(self) -> dict:
        """Return the obs fields to validate against."""
        return self._obs_fields

    def lookup(self, public: bool = False) -> CatLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CatLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, "var_index": self._var_field},
            public=public,
            sources=self._sources,
        )

    def add_new_from(self, key: str, **kwargs):
        """Add validated & new categories.

        Args:
            key: The key referencing the slot in the DataFrame from which to draw terms.
            **kwargs: Additional keyword arguments to pass to create new records
        """
        if key == "all":
            logger.warning(
                "'all' is deprecated, please pass a single key from `.non_validated.keys()` instead!"
            )
            for k in self.non_validated.keys():
                self._cat_vectors[k].add_new(**kwargs)
        else:
            self._cat_vectors[key].add_new(**kwargs)

    @deprecated(new_name="add_new_from('var_index')")
    def add_new_from_var_index(self, **kwargs):
        """Update variable records.

        Args:
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self.add_new_from(key="var_index", **kwargs)

    def validate(self) -> bool:
        """Validate categories.

        This method also registers the validated records in the current instance.

        Returns:
            Whether the AnnData object is validated.
        """
        self._validate_category_error_messages = ""  # reset the error messages

        validated = True
        for _, cat_vector in self._cat_vectors.items():
            cat_vector.validate()
            validated &= cat_vector.is_validated

        self._non_validated = {}  # so it's no longer None
        self._is_validated = validated
        return self._is_validated

    def standardize(self, key: str):
        """Replace synonyms with standardized values.

        Args:
            key: The key referencing the slot in `adata.obs` from which to draw terms. Same as the key in `categoricals`.

                - If "var_index", standardize the var.index.
                - If "all", standardize all obs columns and var.index.

        Inplace modification of the dataset.
        """
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")
        if key == "all":
            logger.warning(
                "'all' is deprecated, please pass a single key from `.non_validated.keys()` instead!"
            )
            for k in self.non_validated.keys():
                self._cat_vectors[k].standardize()
        else:
            self._cat_vectors[key].standardize()


@deprecated(new_name="MuDataCurator")
class MuDataCatManager(CatManager):
    """Categorical manager for `MuData`."""

    def __init__(
        self,
        mdata: MuData | Artifact,
        var_index: dict[str, FieldAttr] | None = None,
        categoricals: dict[str, FieldAttr] | None = None,
        sources: dict[str, SQLRecord] | None = None,
    ) -> None:
        super().__init__(
            dataset=mdata,
            categoricals={},
            sources=sources,
        )
        self._columns_field = (
            var_index or {}
        )  # this is for consistency with BaseCatManager
        self._var_fields = var_index or {}
        self._verify_modality(self._var_fields.keys())
        self._obs_fields = self._parse_categoricals(categoricals or {})
        self._modalities = set(self._var_fields.keys()) | set(self._obs_fields.keys())
        self._obs_df_curator = None
        if "obs" in self._modalities:
            self._obs_df_curator = DataFrameCatManager(
                df=self._dataset.obs,
                columns_field=Feature.name,
                categoricals=self._obs_fields.get("obs", {}),
                sources=self._sources.get("obs"),
            )
        self._mod_adata_curators = {
            modality: AnnDataCatManager(
                data=self._dataset[modality],
                var_index=var_index.get(modality),
                categoricals=self._obs_fields.get(modality),
                sources=self._sources.get(modality),
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
        non_validated = {}
        if (
            self._obs_df_curator is not None
            and len(self._obs_df_curator.non_validated) > 0
        ):
            non_validated["obs"] = self._obs_df_curator.non_validated
        for modality, adata_curator in self._mod_adata_curators.items():
            if len(adata_curator.non_validated) > 0:
                non_validated[modality] = adata_curator.non_validated
        self._non_validated = non_validated
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

    def lookup(self, public: bool = False) -> CatLookup:
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
        return CatLookup(
            categoricals=obs_fields,
            slots={
                **{f"{k}_var_index": v for k, v in self._var_fields.items()},
            },
            public=public,
            sources=self._sources,
        )

    @deprecated(new_name="add_new_from('var_index')")
    def add_new_from_var_index(self, modality: str, **kwargs):
        """Update variable records.

        Args:
            modality: The modality name.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        self._mod_adata_curators[modality].add_new_from(key="var_index", **kwargs)

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
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        modality = modality or "obs"
        if modality in self._mod_adata_curators:
            adata_curator = self._mod_adata_curators[modality]
            adata_curator.add_new_from(key=key, **kwargs)
        if modality == "obs":
            self._obs_df_curator.add_new_from(key=key, **kwargs)
        if key == "var_index":
            self._mod_adata_curators[modality].add_new_from(key=key, **kwargs)

    def validate(self) -> bool:
        """Validate categories."""
        obs_validated = True
        if "obs" in self._modalities:
            logger.info('validating categoricals in "obs"...')
            obs_validated &= self._obs_df_curator.validate()

        mods_validated = True
        for modality, adata_curator in self._mod_adata_curators.items():
            logger.info(f'validating categoricals in modality "{modality}"...')
            mods_validated &= adata_curator.validate()

        self._non_validated = {}  # so it's no longer None
        self._is_validated = obs_validated & mods_validated
        return self._is_validated

    def standardize(self, key: str, modality: str | None = None):
        """Replace synonyms with standardized values.

        Args:
            key: The key referencing the slot in the `MuData`.
            modality: The modality name.

        Inplace modification of the dataset.
        """
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")
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


@deprecated(new_name="SpatialDataCurator")
class SpatialDataCatManager(CatManager):
    """Categorical manager for `SpatialData`."""

    def __init__(
        self,
        sdata: Any,
        var_index: dict[str, FieldAttr],
        categoricals: dict[str, dict[str, FieldAttr]] | None = None,
        sources: dict[str, dict[str, SQLRecord]] | None = None,
        *,
        sample_metadata_key: str | None = "sample",
    ) -> None:
        super().__init__(
            dataset=sdata,
            categoricals={},
            sources=sources,
        )
        if isinstance(sdata, Artifact):
            self._sdata = sdata.load()
        else:
            self._sdata = self._dataset
        self._sample_metadata_key = sample_metadata_key
        self._write_path = None
        self._var_fields = var_index
        self._verify_accessor_exists(self._var_fields.keys())
        self._categoricals = categoricals
        self._table_keys = set(self._var_fields.keys()) | set(
            self._categoricals.keys() - {self._sample_metadata_key}
        )
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

        # check validity of keys in sources
        nonval_keys = []
        for accessor, accessor_sources in self._sources.items():
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
        _maybe_curation_keys_not_present(nonval_keys, "sources")

        # Set up sample level metadata and table Curator objects
        if (
            self._sample_metadata_key is not None
            and self._sample_metadata_key in self._categoricals
        ):
            self._sample_df_curator = DataFrameCatManager(
                df=self._sample_metadata,
                columns_field=Feature.name,
                categoricals=self._categoricals.get(self._sample_metadata_key, {}),
                sources=self._sources.get(self._sample_metadata_key),
            )
        self._table_adata_curators = {
            table: AnnDataCatManager(
                data=self._sdata[table],
                var_index=var_index.get(table),
                categoricals=self._categoricals.get(table),
                sources=self._sources.get(table),
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
    def non_validated(self) -> dict[str, dict[str, list[str]]]:  # type: ignore
        """Return the non-validated features and labels."""
        if self._non_validated is None:
            raise ValidationError("Please run validate() first!")
        non_curated = {}
        if len(self._sample_df_curator.non_validated) > 0:
            non_curated[self._sample_metadata_key] = (
                self._sample_df_curator.non_validated
            )
        for table, adata_curator in self._table_adata_curators.items():
            if len(adata_curator.non_validated) > 0:
                non_curated[table] = adata_curator.non_validated
        return non_curated

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

    def lookup(self, public: bool = False) -> CatLookup:
        """Look up categories.

        Args:
            public: Whether the lookup is performed on the public reference.
        """
        cat_values_dict = list(self.categoricals.values())[0]
        return CatLookup(
            categoricals=cat_values_dict,
            slots={"accessors": cat_values_dict.keys()},
            public=public,
            sources=self._sources,
        )

    @deprecated(new_name="add_new_from('var_index')")
    def add_new_from_var_index(self, table: str, **kwargs) -> None:
        """Save new values from ``.var.index`` of table.

        Args:
            table: The table key.
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        if table in self.non_validated.keys():
            self._table_adata_curators[table].add_new_from(key="var_index", **kwargs)

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
            **kwargs: Additional keyword arguments to pass to create new records.
        """
        if accessor in self.non_validated.keys():
            if accessor in self._table_adata_curators:
                adata_curator = self._table_adata_curators[accessor]
                adata_curator.add_new_from(key=key, **kwargs)
            if accessor == self._sample_metadata_key:
                self._sample_df_curator.add_new_from(key=key, **kwargs)

        if key == "var_index":
            self._table_adata_curators[accessor].add_new_from(key=key, **kwargs)

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
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")

        if accessor == self._sample_metadata_key:
            if key not in self._sample_metadata.columns:
                raise ValueError(f"key '{key}' not present in '{accessor}'!")
        else:
            if (
                key == "var_index" and self._sdata.tables[accessor].var.index is None
            ) or (
                key != "var_index"
                and key not in self._sdata.tables[accessor].obs.columns
            ):
                raise ValueError(f"key '{key}' not present in '{accessor}'!")

        if accessor in self._table_adata_curators.keys():
            adata_curator = self._table_adata_curators[accessor]
            adata_curator.standardize(key)
        if accessor == self._sample_metadata_key:
            self._sample_df_curator.standardize(key)

    def validate(self) -> bool:
        """Validate variables and categorical observations.

        This method also registers the validated records in the current instance:
        - from public sources

        Returns:
            Whether the SpatialData object is validated.
        """
        # add all validated records to the current instance
        sample_validated = True
        if self._sample_df_curator:
            logger.info(f"validating categoricals of '{self._sample_metadata_key}' ...")
            sample_validated &= self._sample_df_curator.validate()

        mods_validated = True
        for table, adata_curator in self._table_adata_curators.items():
            logger.info(f"validating categoricals of table '{table}' ...")
            mods_validated &= adata_curator.validate()

        self._non_validated = {}  # so it's no longer None
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
        """Save the validated SpatialData store and metadata.

        Args:
            description: A description of the dataset.
            key: A path-like key to reference artifact in default storage,
                e.g., `"myartifact.zarr"`. Artifacts with the same key form a version family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
        if not self._is_validated:
            self.validate()
            if not self._is_validated:
                raise ValidationError("Dataset does not validate. Please curate.")

        self._artifact = Artifact.from_spatialdata(
            self._dataset, key=key, description=description, revises=revises, run=run
        ).save()
        return legacy_annotate_artifact(
            self._artifact,
            index_field=self.var_index,
            sample_metadata_key=self._sample_metadata_key,
        )


class TiledbsomaCatManager(CatManager):
    """Categorical manager for `tiledbsoma.Experiment`."""

    def __init__(
        self,
        experiment_uri: UPathStr | Artifact,
        var_index: dict[str, tuple[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        sources: dict[str, SQLRecord] | None = None,
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
        self._sources = sources or {}

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

        # check validity of keys in sources
        valid_arg_keys = valid_obs_keys + valid_var_keys + ["columns"]
        nonval_keys = []
        for arg_key in self._sources.keys():
            if arg_key not in valid_arg_keys:
                nonval_keys.append(arg_key)
        _maybe_curation_keys_not_present(nonval_keys, "sources")

        # register obs columns' names
        register_columns = list(self._obs_fields.keys())
        # register categorical keys as features
        cat_vector = CatVector(
            values_getter=register_columns,
            field=self._columns_field,
            key="columns",
            source=self._sources.get("columns"),
        )
        cat_vector.add_new()

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
                cat_vector = CatVector(
                    values_getter=var_ms_values,
                    field=field,
                    key=var_ms_key,
                    source=self._sources.get(var_ms_key),
                )
                cat_vector.validate()
                non_val = cat_vector._non_validated
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
                cat_vector = CatVector(
                    values_getter=values,
                    field=field,
                    key=key,
                    source=self._sources.get(key),
                )
                cat_vector.validate()
                non_val = cat_vector._non_validated
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
            cat_vector = CatVector(
                values_getter=values,
                field=field,
                key=k,
                source=self._sources.get(k),
            )
            cat_vector.add_new()
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

    def lookup(self, public: bool = False) -> CatLookup:
        """Lookup categories.

        Args:
            public: If "public", the lookup is performed on the public reference.
        """
        return CatLookup(
            categoricals=self._obs_fields,
            slots={"columns": self._columns_field, **self._var_fields_flat},
            public=public,
            sources=self._sources,
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
            cat_vector = CatVector(
                values_getter=values,
                field=field,
                key=k,
                source=self._sources.get(k),
            )
            cat_vector.validate()
            syn_mapper = cat_vector._synonyms
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
                e.g., `"myfolder/mystore.tiledbsoma"`. Artifacts with the same key form a version family.
            revises: Previous version of the artifact. Triggers a revision.
            run: The run that creates the artifact.

        Returns:
            A saved artifact record.
        """
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

        feature_sets = {}
        if len(self._obs_fields) > 0:
            empty_dict = {field.name: [] for field in self._obs_pa_schema}  # type: ignore
            mock_df = pa.Table.from_pydict(
                empty_dict, schema=self._obs_pa_schema
            ).to_pandas()
            # in parallel to https://github.com/laminlabs/lamindb/blob/2a1709990b5736b480c6de49c0ada47fafc8b18d/lamindb/core/_feature_manager.py#L549-L554
            feature_sets["obs"] = Schema.from_df(
                df=mock_df,
                field=self._columns_field,
                mute=True,
            )
        for ms in self._var_fields:
            var_key, var_field = self._var_fields[ms]
            feature_sets[f"{ms}__var"] = Schema.from_values(
                values=self._validated_values[f"{ms}__{var_key}"],
                field=var_field,
                raise_validation_error=False,
            )
        artifact._staged_feature_sets = feature_sets

        feature_ref_is_name = _ref_is_name(self._columns_field)
        features = Feature.lookup().dict()
        for key, field in self._obs_fields.items():
            feature = features.get(key)
            registry = field.field.model
            labels = registry.from_values(
                values=self._validated_values[key],
                field=field,
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


def legacy_annotate_artifact(
    artifact: Artifact,
    *,
    cat_vectors: dict[str, CatVector] | None = None,
    index_field: FieldAttr | dict[str, FieldAttr] | None = None,
    **kwargs,
) -> Artifact:
    from ..models.artifact import add_labels

    if cat_vectors is None:
        cat_vectors = {}

    # annotate with labels
    for key, cat_vector in cat_vectors.items():
        if (
            cat_vector._field.field.model == Feature
            or key == "columns"
            or key == "var_index"
        ):
            continue
        add_labels(
            artifact,
            records=cat_vector.records,
            feature=cat_vector.feature,
            feature_ref_is_name=None,  # do not need anymore
            label_ref_is_name=cat_vector.label_ref_is_name,
            from_curator=True,
        )

    match artifact.otype:
        case "DataFrame":
            artifact.features._add_set_from_df(field=index_field)  # type: ignore
        case "AnnData":
            artifact.features._add_set_from_anndata(  # type: ignore
                var_field=index_field,
            )
        case "MuData":
            artifact.features._add_set_from_mudata(var_fields=index_field)  # type: ignore
        case "SpatialData":
            artifact.features._add_set_from_spatialdata(  # type: ignore
                sample_metadata_key=kwargs.get("sample_metadata_key", "sample"),
                var_fields=index_field,
            )
        case _:
            raise NotImplementedError  # pragma: no cover

    return artifact


# backward compat constructors ------------------


@classmethod  # type: ignore
def from_df(
    cls,
    df: pd.DataFrame,
    categoricals: dict[str, FieldAttr] | None = None,
    columns: FieldAttr = Feature.name,
    organism: str | None = None,
) -> DataFrameCatManager:
    if organism is not None:
        logger.warning("organism is ignored, define it on the dtype level")
    return DataFrameCatManager(
        df=df,
        categoricals=categoricals,
        columns_field=columns,
    )


@classmethod  # type: ignore
def from_anndata(
    cls,
    data: ad.AnnData | UPathStr,
    var_index: FieldAttr,
    categoricals: dict[str, FieldAttr] | None = None,
    obs_columns: FieldAttr = Feature.name,
    organism: str | None = None,
    sources: dict[str, SQLRecord] | None = None,
) -> AnnDataCatManager:
    if organism is not None:
        logger.warning("organism is ignored, define it on the dtype level")
    return AnnDataCatManager(
        data=data,
        var_index=var_index,
        categoricals=categoricals,
        obs_columns=obs_columns,
        sources=sources,
    )


@classmethod  # type: ignore
def from_mudata(
    cls,
    mdata: MuData | UPathStr,
    var_index: dict[str, dict[str, FieldAttr]],
    categoricals: dict[str, FieldAttr] | None = None,
    organism: str | None = None,
) -> MuDataCatManager:
    if not is_package_installed("mudata"):
        raise ImportError("Please install mudata: pip install mudata")
    if organism is not None:
        logger.warning("organism is ignored, define it on the dtype level")
    return MuDataCatManager(
        mdata=mdata,
        var_index=var_index,
        categoricals=categoricals,
    )


@classmethod  # type: ignore
def from_tiledbsoma(
    cls,
    experiment_uri: UPathStr,
    var_index: dict[str, tuple[str, FieldAttr]],
    categoricals: dict[str, FieldAttr] | None = None,
    obs_columns: FieldAttr = Feature.name,
    organism: str | None = None,
    sources: dict[str, SQLRecord] | None = None,
) -> TiledbsomaCatManager:
    if organism is not None:
        logger.warning("organism is ignored, define it on the dtype level")
    return TiledbsomaCatManager(
        experiment_uri=experiment_uri,
        var_index=var_index,
        categoricals=categoricals,
        obs_columns=obs_columns,
        sources=sources,
    )


@classmethod  # type: ignore
def from_spatialdata(
    cls,
    sdata: SpatialData | UPathStr,
    var_index: dict[str, FieldAttr],
    categoricals: dict[str, dict[str, FieldAttr]] | None = None,
    organism: str | None = None,
    sources: dict[str, dict[str, SQLRecord]] | None = None,
    *,
    sample_metadata_key: str = "sample",
):
    if not is_package_installed("spatialdata"):
        raise ImportError("Please install spatialdata: pip install spatialdata")
    if organism is not None:
        logger.warning("organism is ignored, define it on the dtype level")
    return SpatialDataCatManager(
        sdata=sdata,
        var_index=var_index,
        categoricals=categoricals,
        sources=sources,
        sample_metadata_key=sample_metadata_key,
    )


CatManager.from_df = from_df  # type: ignore
CatManager.from_anndata = from_anndata  # type: ignore
CatManager.from_mudata = from_mudata  # type: ignore
CatManager.from_spatialdata = from_spatialdata  # type: ignore
CatManager.from_tiledbsoma = from_tiledbsoma  # type: ignore
