"""Curators.

.. versionadded:: 1.1.0

.. autosummary::
   :toctree: .

   Curator
   DataFrameCurator
   SlotsCurator
   AnnDataCurator
   MuDataCurator
   SpatialDataCurator

CatManager:

.. autosummary::
   :toctree: .

   CatManager
   DataFrameCatManager
   AnnDataCatManager
   MuDataCatManager
   TiledbsomaCatManager
   SpatialDataCatManager
   CurateLookup

"""

from __future__ import annotations

import copy
import re
from itertools import chain
from typing import TYPE_CHECKING, Any, Literal

import anndata as ad
import lamindb_setup as ln_setup
import pandas as pd
import pandera
import pyarrow as pa
from lamin_utils import colors, logger
from lamindb_setup.core import deprecated, upath
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.upath import UPath

from lamindb.core.storage._backed_access import backed_access

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr
    from mudata import MuData
    from spatialdata import SpatialData

    from lamindb.core.types import ScverseDataStructures
    from lamindb.models import Record
from lamindb.base.types import FieldAttr  # noqa
from lamindb.core._settings import settings
from lamindb.models import (
    Artifact,
    Feature,
    Record,
    Run,
    Schema,
    ULabel,
)
from lamindb.models.artifact import (
    add_labels,
    data_is_anndata,
    data_is_mudata,
    data_is_spatialdata,
)
from lamindb.models.feature import parse_dtype, parse_dtype_single_cat
from lamindb.models._from_values import _format_values

from ..errors import InvalidArgument, ValidationError
from anndata import AnnData

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableMapping
    from typing import Any

    from lamindb_setup.core.types import UPathStr

    from lamindb.models.query_set import RecordList


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

    Example::

        curator = ln.curators.DataFrameCurator(...)
        curator.cat.lookup()["cell_type"].alveolar_type_1_fibroblast_cell

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
        else:  # pragma: no cover
            return colors.warning("No fields are found!")


CAT_MANAGER_DOCSTRING = """Manage categoricals by updating registries."""


SLOTS_DOCSTRING = """Curator objects by slot.

.. versionadded:: 1.1.1
"""


VALIDATE_DOCSTRING = """Validate dataset against Schema.

Raises:
    lamindb.errors.ValidationError: If validation fails.
"""

SAVE_ARTIFACT_DOCSTRING = """Save an annotated artifact.

Args:
    key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a version family.
    description: A description.
    revises: Previous version of the artifact. Is an alternative way to passing `key` to trigger a new version.
    run: The run that creates the artifact.

Returns:
    A saved artifact record.
"""


class Curator:
    """Dataset curator.

    A `Curator` object makes it easy to validate, standardize & annotate datasets.

    .. versionadded:: 1.1.0

    See:
        - :class:`~lamindb.curators.DataFrameCurator`
        - :class:`~lamindb.curators.AnnDataCurator`
        - :class:`~lamindb.curators.MuDataCurator`
        - :class:`~lamindb.curators.SpatialDataCurator`
    """

    def __init__(self, dataset: Any, schema: Schema | None = None):
        self._artifact: Artifact = None  # pass the dataset as an artifact
        self._dataset: Any = dataset  # pass the dataset as a UPathStr or data object
        if isinstance(self._dataset, Artifact):
            self._artifact = self._dataset
            if self._artifact.otype in {
                "DataFrame",
                "AnnData",
                "MuData",
                "SpatialData",
            }:
                self._dataset = self._dataset.load()
        self._schema: Schema | None = schema
        self._is_validated: bool = False
        self._cat_manager: CatManager = None  # is None for CatManager curators

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> bool | str:
        """{}"""  # noqa: D415
        pass  # pragma: no cover

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
        # Note that this docstring has to be consistent with the Artifact()
        # constructor signature
        pass


class SlotsCurator(Curator):
    """Curator for a dataset with slots.

    Args:
        dataset: The dataset to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    .. versionadded:: 1.3.0
    """

    def __init__(
        self,
        dataset: Any,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        self._slots: dict[str, DataFrameCurator] = {}

        # used in MuDataCurator and SpatialDataCurator
        # in form of {table/modality_key: var_field}
        self._var_fields: dict[str, FieldAttr] = {}
        # in form of {table/modality_key: categoricals}
        self._categoricals: dict[str, dict[str, FieldAttr]] = {}

    @property
    @doc_args(SLOTS_DOCSTRING)
    def slots(self) -> dict[str, DataFrameCurator]:
        """{}"""  # noqa: D415
        return self._slots

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}"""  # noqa: D415
        for _, curator in self._slots.items():
            curator.validate()

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
        if not self._is_validated:
            self.validate()

        # default implementation for MuDataCurator and SpatialDataCurator
        return save_artifact(  # type: ignore
            self._dataset,
            key=key,
            description=description,
            fields=self._categoricals,
            index_field=self._var_fields,
            artifact=self._artifact,
            revises=revises,
            run=run,
            schema=self._schema,
        )


class DataFrameCurator(Curator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for a DataFrame object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    .. versionadded:: 1.1.0

    Args:
        dataset: The DataFrame-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Example::

        import lamindb as ln
        import bionty as bt

        # define valid labels
        perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
        ln.ULabel(name="DMSO", type=perturbation).save()
        ln.ULabel(name="IFNG", type=perturbation).save()
        bt.CellType.from_source(name="B cell").save()
        bt.CellType.from_source(name="T cell").save()

        # define schema
        schema = ln.Schema(
            name="small_dataset1_obs_level_metadata",
            features=[
                ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
                ln.Feature(name="sample_note", dtype=str).save(),
                ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
                ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
            ],
        ).save()

        # curate a DataFrame
        df = datasets.small_dataset1(otype="DataFrame")
        curator = ln.curators.DataFrameCurator(df, schema)
        artifact = curator.save_artifact(key="example_datasets/dataset1.parquet")
        assert artifact.schema == schema
    """

    def __init__(
        self,
        dataset: pd.DataFrame | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        categoricals = {}
        if schema.n > 0:
            # populate features
            pandera_columns = {}
            for feature in schema.features.all():
                pandera_dtype = (
                    feature.dtype if not feature.dtype.startswith("cat") else "category"
                )
                pandera_columns[feature.name] = pandera.Column(
                    pandera_dtype, nullable=feature.nullable
                )
                if feature.dtype.startswith("cat"):
                    categoricals[feature.name] = parse_dtype(feature.dtype)[0]["field"]
            self._pandera_schema = pandera.DataFrameSchema(
                pandera_columns, coerce=schema.coerce_dtype
            )
        else:
            assert schema.itype is not None  # noqa: S101
        self._cat_manager = DataFrameCatManager(
            self._dataset,
            columns=parse_dtype_single_cat(schema.itype, is_itype=True)["field"],
            categoricals=categoricals,
        )

    @property
    @doc_args(CAT_MANAGER_DOCSTRING)
    def cat(self) -> CatManager:
        """{}"""  # noqa: D415
        return self._cat_manager

    def standardize(self) -> None:
        """Standardize the dataset.

        - Adds missing columns for features
        - Fills missing values for features with default values
        """
        for feature in self._schema.members:
            if feature.name not in self._dataset.columns:
                if feature.default_value is not None or feature.nullable:
                    fill_value = (
                        feature.default_value
                        if feature.default_value is not None
                        else pd.NA
                    )
                    if feature.dtype.startswith("cat"):
                        self._dataset[feature.name] = pd.Categorical(
                            [fill_value] * len(self._dataset)
                        )
                    else:
                        self._dataset[feature.name] = fill_value
                    logger.important(
                        f"added column {feature.name} with fill value {fill_value}"
                    )
                else:
                    raise ValidationError(
                        f"Missing column {feature.name} cannot be added because is not nullable and has no default value"
                    )
            else:
                if feature.default_value is not None:
                    if isinstance(
                        self._dataset[feature.name].dtype, pd.CategoricalDtype
                    ):
                        if (
                            feature.default_value
                            not in self._dataset[feature.name].cat.categories
                        ):
                            self._dataset[feature.name] = self._dataset[
                                feature.name
                            ].cat.add_categories(feature.default_value)
                    self._dataset[feature.name] = self._dataset[feature.name].fillna(
                        feature.default_value
                    )

    def _cat_manager_validate(self) -> None:
        self._cat_manager.validate()
        if self._cat_manager._is_validated:
            self._is_validated = True
        else:
            self._is_validated = False
            raise ValidationError(self._cat_manager._validate_category_error_messages)

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}"""  # noqa: D415
        if self._schema.n > 0:
            try:
                # first validate through pandera
                self._pandera_schema.validate(self._dataset)
                # then validate lamindb categoricals
                self._cat_manager_validate()
            except pandera.errors.SchemaError as err:
                self._is_validated = False
                # .exconly() doesn't exist on SchemaError
                raise ValidationError(str(err)) from err
        else:
            self._cat_manager_validate()

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
        if not self._is_validated:
            self.validate()  # raises ValidationError if doesn't validate
        result = parse_dtype_single_cat(self._schema.itype, is_itype=True)
        return save_artifact(  # type: ignore
            self._dataset,
            description=description,
            fields=self._cat_manager.categoricals,
            index_field=result["field"],
            key=key,
            artifact=self._artifact,
            revises=revises,
            run=run,
            schema=self._schema,
        )


class AnnDataCurator(SlotsCurator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for an AnnData object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    .. versionadded:: 1.1.0

    Args:
        dataset: The AnnData-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Example::

        import lamindb as ln
        import bionty as bt

        # define valid labels
        perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
        ln.ULabel(name="DMSO", type=perturbation).save()
        ln.ULabel(name="IFNG", type=perturbation).save()
        bt.CellType.from_source(name="B cell").save()
        bt.CellType.from_source(name="T cell").save()

        # define obs schema
        obs_schema = ln.Schema(
            name="small_dataset1_obs_level_metadata",
            features=[
                ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
                ln.Feature(name="sample_note", dtype=str).save(),
                ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
                ln.Feature(name="cell_type_by_model", dtype=bt.CellType").save(),
            ],
        ).save()

        # define var schema
        var_schema = ln.Schema(
            name="scRNA_seq_var_schema",
            itype=bt.Gene.ensembl_gene_id,
            dtype=int,
        ).save()

        # define composite schema
        anndata_schema = ln.Schema(
            name="small_dataset1_anndata_schema",
            otype="AnnData",
            components={"obs": obs_schema, "var": var_schema},
        ).save()

        # curate an AnnData
        adata = ln.core.datasets.small_dataset1(otype="AnnData")
        curator = ln.curators.AnnDataCurator(adata, anndata_schema)
        artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
        assert artifact.schema == anndata_schema
    """

    def __init__(
        self,
        dataset: AnnData | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_anndata(self._dataset):
            raise InvalidArgument("dataset must be AnnData-like.")
        if schema.otype != "AnnData":
            raise InvalidArgument("Schema otype must be 'AnnData'.")
        # TODO: also support slots other than obs and var
        self._slots = {
            slot: DataFrameCurator(
                (
                    getattr(self._dataset, slot).T
                    if slot == "var"
                    else getattr(self._dataset, slot)
                ),
                slot_schema,
            )
            for slot, slot_schema in schema.slots.items()
            if slot in {"obs", "var"}
        }

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
        if not self._is_validated:
            self.validate()
        return save_artifact(  # type: ignore
            self._dataset,
            description=description,
            fields=self.slots["obs"]._cat_manager.categoricals,
            index_field=(
                parse_dtype_single_cat(self.slots["var"]._schema.itype, is_itype=True)[
                    "field"
                ]
                if "var" in self._slots
                else None
            ),
            key=key,
            artifact=self._artifact,
            revises=revises,
            run=run,
            schema=self._schema,
        )


def _assign_var_fields_categoricals_multimodal(
    modality: str | None,
    slot_type: str,
    slot: str,
    slot_schema: Schema,
    var_fields: dict[str, FieldAttr],
    categoricals: dict[str, dict[str, FieldAttr]],
    slots: dict[str, DataFrameCurator],
) -> None:
    """Assigns var_fields and categoricals for multimodal data curators."""
    if modality is not None:
        # Makes sure that all tables are present
        var_fields[modality] = None
        categoricals[modality] = {}

    if slot_type == "var":
        var_field = parse_dtype_single_cat(slot_schema.itype, is_itype=True)["field"]
        if modality is None:
            # This should rarely/never be used since tables should have different var fields
            var_fields[slot] = var_field  # pragma: no cover
        else:
            # Note that this is NOT nested since the nested key is always "var"
            var_fields[modality] = var_field
    else:
        obs_fields = slots[slot]._cat_manager.categoricals
        if modality is None:
            categoricals[slot] = obs_fields
        else:
            # Note that this is NOT nested since the nested key is always "obs"
            categoricals[modality] = obs_fields


class MuDataCurator(SlotsCurator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for a MuData object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    .. versionadded:: 1.3.0

    Args:
        dataset: The MuData-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Example::

        import lamindb as ln
        import bionty as bt

        # define the global obs schema
        obs_schema = ln.Schema(
            name="mudata_papalexi21_subset_obs_schema",
            features=[
                ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
                ln.Feature(name="replicate", dtype="cat[ULabel[Replicate]]").save(),
            ],
        ).save()

        # define the ['rna'].obs schema
        obs_schema_rna = ln.Schema(
            name="mudata_papalexi21_subset_rna_obs_schema",
            features=[
                ln.Feature(name="nCount_RNA", dtype=int).save(),
                ln.Feature(name="nFeature_RNA", dtype=int).save(),
                ln.Feature(name="percent.mito", dtype=float).save(),
            ],
            coerce_dtype=True,
        ).save()

        # define the ['hto'].obs schema
        obs_schema_hto = ln.Schema(
            name="mudata_papalexi21_subset_hto_obs_schema",
            features=[
                ln.Feature(name="nCount_HTO", dtype=int).save(),
                ln.Feature(name="nFeature_HTO", dtype=int).save(),
                ln.Feature(name="technique", dtype=bt.ExperimentalFactor).save(),
            ],
            coerce_dtype=True,
        ).save()

        # define ['rna'].var schema
        var_schema_rna = ln.Schema(
            name="mudata_papalexi21_subset_rna_var_schema",
            itype=bt.Gene.symbol,
            dtype=float,
        ).save()

        # define composite schema
        mudata_schema = ln.Schema(
            name="mudata_papalexi21_subset_mudata_schema",
            otype="MuData",
            components={
                "obs": obs_schema,
                "rna:obs": obs_schema_rna,
                "hto:obs": obs_schema_hto,
                "rna:var": var_schema_rna,
            },
        ).save()

        # curate a MuData
        mdata = ln.core.datasets.mudata_papalexi21_subset()
        bt.settings.organism = "human" # set the organism
        curator = ln.curators.MuDataCurator(mdata, mudata_schema)
        artifact = curator.save_artifact(key="example_datasets/mudata_papalexi21_subset.h5mu")
        assert artifact.schema == mudata_schema
    """

    def __init__(
        self,
        dataset: MuData | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_mudata(self._dataset):
            raise InvalidArgument("dataset must be MuData-like.")
        if schema.otype != "MuData":
            raise InvalidArgument("Schema otype must be 'MuData'.")

        for slot, slot_schema in schema.slots.items():
            # Assign to _slots
            if ":" in slot:
                modality, modality_slot = slot.split(":")
                schema_dataset = self._dataset.__getitem__(modality)
            else:
                modality, modality_slot = None, slot
                schema_dataset = self._dataset
            self._slots[slot] = DataFrameCurator(
                (
                    getattr(schema_dataset, modality_slot).T
                    if modality_slot == "var"
                    else getattr(schema_dataset, modality_slot)
                ),
                slot_schema,
            )
            _assign_var_fields_categoricals_multimodal(
                modality=modality,
                slot_type=modality_slot,
                slot=slot,
                slot_schema=slot_schema,
                var_fields=self._var_fields,
                categoricals=self._categoricals,
                slots=self._slots,
            )

        # for consistency with BaseCatManager
        self._columns_field = self._var_fields


class SpatialDataCurator(SlotsCurator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for a SpatialData object.

    See also :class:`~lamindb.Curator` and :class:`~lamindb.Schema`.

    .. versionadded:: 1.3.0

    Args:
        dataset: The SpatialData-like object to validate & annotate.
        schema: A `Schema` object that defines the validation constraints.

    Example::

        import lamindb as ln
        import bionty as bt

        # define sample schema
        sample_schema = ln.Schema(
            name="blobs_sample_level_metadata",
            features=[
                ln.Feature(name="assay", dtype=bt.ExperimentalFactor).save(),
                ln.Feature(name="disease", dtype=bt.Disease).save(),
                ln.Feature(name="development_stage", dtype=bt.DevelopmentalStage).save(),
            ],
            coerce_dtype=True
        ).save()

        # define table obs schema
        blobs_obs_schema = ln.Schema(
            name="blobs_obs_level_metadata",
            features=[
                ln.Feature(name="sample_region", dtype="str").save(),
            ],
            coerce_dtype=True
        ).save()

        # define table var schema
        blobs_var_schema = ln.Schema(
            name="blobs_var_schema",
            itype=bt.Gene.ensembl_gene_id,
            dtype=int
        ).save()

        # define composite schema
        spatialdata_schema = ln.Schema(
            name="blobs_spatialdata_schema",
            otype="SpatialData",
            components={
                "sample": sample_schema,
                "table:obs": blobs_obs_schema,
                "table:var": blobs_var_schema,
        }).save()

        # curate a SpatialData
        spatialdata = ln.core.datasets.spatialdata_blobs()
        curator = ln.curators.SpatialDataCurator(spatialdata, spatialdata_schema)
        try:
            curator.validate()
        except ln.errors.ValidationError as error:
            print(error)

        # validate again (must pass now) and save artifact
        artifact = curator.save_artifact(key="example_datasets/spatialdata1.zarr")
        assert artifact.schema == spatialdata_schema
    """

    def __init__(
        self,
        dataset: SpatialData | Artifact,
        schema: Schema,
        *,
        sample_metadata_key: str | None = "sample",
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_spatialdata(self._dataset):
            raise InvalidArgument("dataset must be SpatialData-like.")
        if schema.otype != "SpatialData":
            raise InvalidArgument("Schema otype must be 'SpatialData'.")

        for slot, slot_schema in schema.slots.items():
            # Assign to _slots
            if ":" in slot:
                table_key, table_slot = slot.split(":")
                schema_dataset = self._dataset.tables.__getitem__(table_key)
            # sample metadata (does not have a `:` separator)
            else:
                table_key = None
                table_slot = slot
                schema_dataset = self._dataset.get_attrs(
                    key=sample_metadata_key, return_as="df", flatten=True
                )

            self._slots[slot] = DataFrameCurator(
                (
                    getattr(schema_dataset, table_slot).T
                    if table_slot == "var"
                    else (
                        getattr(schema_dataset, table_slot)
                        if table_slot != sample_metadata_key
                        else schema_dataset
                    )  # just take the schema_dataset if it's the sample metadata key
                ),
                slot_schema,
            )

            _assign_var_fields_categoricals_multimodal(
                modality=table_key,
                slot_type=table_slot,
                slot=slot,
                slot_schema=slot_schema,
                var_fields=self._var_fields,
                categoricals=self._categoricals,
                slots=self._slots,
            )

        # for consistency with BaseCatManager
        self._columns_field = self._var_fields


class CatManager:
    """Manage valid categoricals by updating registries.

    A `CatManager` object makes it easy to validate, standardize & annotate datasets.

    Example::

        cat_manager = ln.CatManager(
            dataset,
            # define validation criteria as mappings
            columns=Feature.name,  # map column names
            categoricals={"perturbation": ULabel.name},  # map categories
        )
        cat_manager.validate()  # validate the dataframe
        artifact = cat_manager.save_artifact(description="my RNA-seq")
        artifact.describe()  # see annotations

    `cat_manager.validate()` maps values within `df` according to the mapping criteria and logs validated & problematic values.

    If you find non-validated values, you have several options:

    - new values found in the data can be registered using :meth:`~lamindb.curators.DataFrameCatManager.add_new_from`
    - non-validated values can be accessed using :meth:`~lamindb.curators.DataFrameCatManager.non_validated` and addressed manually
    """

    def __init__(self, *, dataset, categoricals, sources, organism, columns_field=None):
        # the below is shared with Curator
        self._artifact: Artifact = None  # pass the dataset as an artifact
        self._dataset: Any = dataset  # pass the dataset as a UPathStr or data object
        if isinstance(self._dataset, Artifact):
            self._artifact = self._dataset
            if self._artifact.otype in {"DataFrame", "AnnData"}:
                self._dataset = self._dataset.load()
        self._is_validated: bool = False
        # shared until here
        self._categoricals = categoricals or {}
        self._non_validated = None
        self._organism = organism
        self._sources = sources or {}
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

    def validate(self) -> bool:
        """Validate dataset.

        This method also registers the validated records in the current instance.

        Returns:
            The boolean `True` if the dataset is validated. Otherwise, a string with the error message.
        """
        pass

    def standardize(self, key: str) -> None:
        """Replace synonyms with standardized values.

        Inplace modification of the dataset.

        Args:
            key: The name of the column to standardize.

        Returns:
            None
        """
        pass  # pragma: no cover

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self,
        *,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ) -> Artifact:
        """{}"""  # noqa: D415
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
                key=key,
                description=description,
                fields=self.categoricals,
                index_field=self._columns_field,
                artifact=self._artifact,
                revises=revises,
                run=run,
                schema=None,
                organism=self._organism,
            )
        finally:
            settings.verbosity = verbosity

        return self._artifact


class DataFrameCatManager(CatManager):
    """Curation flow for a DataFrame object.

    See also :class:`~lamindb.Curator`.

    Args:
        df: The DataFrame object to curate.
        columns: The field attribute for the feature column.
        categoricals: A dictionary mapping column names to registry_field.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping column names to Source records.

    Returns:
        A curator object.

    Example::

        import bionty as bt
        curator = ln.Curator.from_df(
            df,
            categoricals={
                "cell_type_ontology_id": bt.CellType.ontology_id,
                "donor_id": ULabel.name
            }
        )
    """

    def __init__(
        self,
        df: pd.DataFrame | Artifact,
        columns: FieldAttr = Feature.name,
        categoricals: dict[str, FieldAttr] | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
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
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")
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


class AnnDataCatManager(CatManager):
    """Manage categorical curation.

    Args:
        data: The AnnData object or an AnnData-like path.
        var_index: The registry field for mapping the ``.var`` index.
        categoricals: A dictionary mapping ``.obs.columns`` to a registry field.
        obs_columns: The registry field for mapping the ``.obs.columns``.
        verbosity: The verbosity level.
        organism: The organism name.
        sources: A dictionary mapping ``.obs.columns`` to Source records.

    Example::

        import bionty as bt
        curator = ln.Curator.from_anndata(
            adata,
            var_index=bt.Gene.ensembl_gene_id,
            categoricals={
                "cell_type_ontology_id": bt.CellType.ontology_id,
                "donor_id": ULabel.name
            },
            organism="human",
        )
    """

    def __init__(
        self,
        data: ad.AnnData | Artifact,
        var_index: FieldAttr | None = None,
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> None:
        if isinstance(var_index, str):
            raise TypeError("var_index parameter has to be a bionty field")

        if not data_is_anndata(data):
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
            organism=organism,
            columns_field=var_index,
        )
        self._adata = self._dataset
        self._obs_df_curator = DataFrameCatManager(
            df=self._adata.obs,
            categoricals=self.categoricals,
            columns=obs_columns,
            verbosity=verbosity,
            organism=None,
            sources=self._sources,
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
        if self.var_index is not None:
            update_registry(
                values=list(self._adata.var.index),
                field=self.var_index,
                key="var_index",
                validated_only=validated_only,
                organism=self._organism,
                source=self._sources.get("var_index"),
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
        if self.var_index is not None:
            validated_var, non_validated_var = validate_categories(
                self._adata.var.index,
                field=self._var_field,
                key="var_index",
                source=self._sources.get("var_index"),
                hint_print=".add_new_from_var_index()",
                organism=self._organism,  # type: ignore
            )
        else:
            validated_var = True
            non_validated_var = []
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
        if self._artifact is not None:
            raise RuntimeError("can't mutate the dataset when an artifact is passed!")
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


class MuDataCatManager(CatManager):
    """Curation flow for a ``MuData`` object.

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

    Example::

        import bionty as bt
        curator = ln.Curator.from_mudata(
            mdata,
            var_index={
                "rna": bt.Gene.ensembl_gene_id,
                "adt": CellMarker.name
            },
            categoricals={
                "cell_type_ontology_id": bt.CellType.ontology_id,
                "donor_id": ULabel.name
            },
            organism="human",
        )
    """

    def __init__(
        self,
        mdata: MuData | Artifact,
        var_index: dict[str, FieldAttr] | None = None,
        categoricals: dict[str, FieldAttr] | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
    ) -> None:
        super().__init__(
            dataset=mdata,
            categoricals={},
            sources=sources,
            organism=organism,
        )
        self._columns_field = (
            var_index or {}
        )  # this is for consistency with BaseCatManager
        self._var_fields = var_index or {}
        self._verify_modality(self._var_fields.keys())
        self._obs_fields = self._parse_categoricals(categoricals or {})
        self._modalities = set(self._var_fields.keys()) | set(self._obs_fields.keys())
        self._verbosity = verbosity
        self._obs_df_curator = None
        if "obs" in self._modalities:
            self._obs_df_curator = DataFrameCatManager(
                df=self._dataset.obs,
                columns=Feature.name,
                categoricals=self._obs_fields.get("obs", {}),
                verbosity=verbosity,
                sources=self._sources.get("obs"),
                organism=organism,
            )
        self._mod_adata_curators = {
            modality: AnnDataCatManager(
                data=self._dataset[modality],
                var_index=var_index.get(modality),
                categoricals=self._obs_fields.get(modality),
                verbosity=verbosity,
                sources=self._sources.get(modality),
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


class SpatialDataCatManager(CatManager):
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
        verbosity: The verbosity level of the logger.
        sample_metadata_key: The key in ``.attrs`` that stores the sample level metadata.

    Example::

        import bionty as bt
        curator = SpatialDataCatManager(
            sdata,
            var_index={
                "table_1": bt.Gene.ensembl_gene_id,
            },
            categoricals={
                "table1":
                    {"cell_type_ontology_id": bt.CellType.ontology_id, "donor_id": ULabel.name},
                "sample":
                    {"experimental_factor": bt.ExperimentalFactor.name},
            },
            organism="human",
        )
    """

    def __init__(
        self,
        sdata: Any,
        var_index: dict[str, FieldAttr],
        categoricals: dict[str, dict[str, FieldAttr]] | None = None,
        verbosity: str = "hint",
        organism: str | None = None,
        sources: dict[str, dict[str, Record]] | None = None,
        *,
        sample_metadata_key: str | None = "sample",
    ) -> None:
        super().__init__(
            dataset=sdata,
            categoricals={},
            sources=sources,
            organism=organism,
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
                columns=Feature.name,
                categoricals=self._categoricals.get(self._sample_metadata_key, {}),
                verbosity=verbosity,
                sources=self._sources.get(self._sample_metadata_key),
                organism=organism,
            )
        self._table_adata_curators = {
            table: AnnDataCatManager(
                data=self._sdata[table],
                var_index=var_index.get(table),
                categoricals=self._categoricals.get(table),
                verbosity=verbosity,
                sources=self._sources.get(table),
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
    def non_validated(self) -> dict[str, dict[str, list[str]]]:  # type: ignore
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
            if "var_index" in self._non_validated[table]:
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
                f"Accessor {accessor} is not in 'categoricals'. Include it when creating the SpatialDataCatManager."
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

        return save_artifact(
            self._sdata,
            description=description,
            fields=self.categoricals,
            index_field=self.var_index,
            key=key,
            artifact=self._artifact,
            revises=revises,
            run=run,
            schema=None,
            organism=self._organism,
            sample_metadata_key=self._sample_metadata_key,
        )


class TiledbsomaCatManager(CatManager):
    """Curation flow for `tiledbsoma.Experiment`.

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

    Example::

        import bionty as bt
        curator = ln.Curator.from_tiledbsoma(
            "./my_array_store.tiledbsoma",
            var_index={"RNA": ("var_id", bt.Gene.symbol)},
            categoricals={
                "cell_type_ontology_id": bt.CellType.ontology_id,
                "donor_id": ULabel.name
            },
            organism="human",
        )
    """

    def __init__(
        self,
        experiment_uri: UPathStr | Artifact,
        var_index: dict[str, tuple[str, FieldAttr]],
        categoricals: dict[str, FieldAttr] | None = None,
        obs_columns: FieldAttr = Feature.name,
        organism: str | None = None,
        sources: dict[str, Record] | None = None,
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
                )
                _, non_val = validate_categories(
                    values=var_ms_values,
                    field=field,
                    key=var_ms_key,
                    organism=organism,
                    source=self._sources.get(var_ms_key),
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
                )
                _, non_val = validate_categories(
                    values=values,
                    field=field,
                    key=key,
                    organism=organism,
                    source=self._sources.get(key),
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
            organism = check_registry_organism(
                self._columns_field.field.model, self._organism
            ).get("organism")
            empty_dict = {field.name: [] for field in self._obs_pa_schema}  # type: ignore
            mock_df = pa.Table.from_pydict(
                empty_dict, schema=self._obs_pa_schema
            ).to_pandas()
            # in parallel to https://github.com/laminlabs/lamindb/blob/2a1709990b5736b480c6de49c0ada47fafc8b18d/lamindb/core/_feature_manager.py#L549-L554
            feature_sets["obs"] = Schema.from_df(
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
            feature_sets[f"{ms}__var"] = Schema.from_values(
                values=self._validated_values[f"{ms}__{var_key}"],
                field=var_field,
                organism=organism,
                raise_validation_error=False,
            )
        artifact._staged_feature_sets = feature_sets

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


def _restrict_obs_fields(
    obs: pd.DataFrame, obs_fields: dict[str, FieldAttr]
) -> dict[str, FieldAttr]:
    """Restrict the obs fields only available obs fields.

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


class CellxGeneAnnDataCatManager(AnnDataCatManager):
    """Annotation flow of AnnData based on CELLxGENE schema."""

    categoricals_defaults = {
        "cell_type": "unknown",
        "development_stage": "unknown",
        "disease": "normal",
        "donor_id": "unknown",
        "self_reported_ethnicity": "unknown",
        "sex": "unknown",
        "suspension_type": "cell",
        "tissue_type": "tissue",
    }

    def __init__(
        self,
        adata: ad.AnnData,
        categoricals: dict[str, FieldAttr] | None = None,
        organism: Literal["human", "mouse"] = "human",
        *,
        defaults: dict[str, str] = None,
        extra_sources: dict[str, Record] = None,
        schema_version: Literal["4.0.0", "5.0.0", "5.1.0", "5.2.0"] = "5.2.0",
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
            schema_version: The CELLxGENE schema version to curate against.
            verbosity: The verbosity level.

        """
        import bionty as bt

        from ._cellxgene_schemas import (
            _create_sources,
            _init_categoricals_additional_values,
        )

        # Add defaults first to ensure that we fetch valid sources
        if defaults:
            _add_defaults_to_obs(adata.obs, defaults)

        # Filter categoricals based on what's present in adata
        if categoricals is None:
            categoricals = self._get_cxg_categoricals()
        categoricals = _restrict_obs_fields(adata.obs, categoricals)

        # Configure sources
        self.sources = _create_sources(categoricals, schema_version, organism)
        self.schema_version = schema_version
        self.schema_reference = f"https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/{schema_version}/schema.md"
        # These sources are not a part of the cellxgene schema but rather passed through.
        # This is useful when other Curators extend the CELLxGENE curator
        if extra_sources:
            self.sources = self.sources | extra_sources

        _init_categoricals_additional_values()

        super().__init__(
            data=adata,
            var_index=bt.Gene.ensembl_gene_id,
            categoricals=categoricals,
            verbosity=verbosity,
            organism=organism,
            sources=self.sources,
        )

    @classmethod
    @deprecated(new_name="categoricals_defaults")
    def _get_categoricals_defaults(cls) -> dict[str, str]:
        return cls.categoricals_defaults

    @classmethod
    def _get_cxg_categoricals(cls) -> dict[str, FieldAttr]:
        """Returns the CELLxGENE required fields."""
        from ._cellxgene_schemas import _get_cxg_categoricals

        return _get_cxg_categoricals()

    def validate(self) -> bool:
        """Validates the AnnData object against most cellxgene requirements."""
        from ._cellxgene_schemas import RESERVED_NAMES

        # Verify that all required obs columns are present
        required_columns = list(self.categoricals_defaults.keys()) + ["donor_id"]
        missing_obs_fields = [
            name
            for name in required_columns
            if name not in self._adata.obs.columns
            and f"{name}_ontology_term_id" not in self._adata.obs.columns
        ]
        if len(missing_obs_fields) > 0:
            logger.error(
                f"missing required obs columns {_format_values(missing_obs_fields)}\n"
                "    → consider initializing a Curate object with `defaults=cxg.CellxGeneAnnDataCatManager.categoricals_defaults` to automatically add these columns with default values"
            )
            return False

        # Verify that no cellxgene reserved names are present
        matched_columns = [
            column for column in self._adata.obs.columns if column in RESERVED_NAMES
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

        def _convert_name_to_ontology_id(values: pd.Series, field: FieldAttr):
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
                mapped_column = _convert_name_to_ontology_id(
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
    def validate_values(cls, values: pd.Series) -> list[str]:
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
    def validate_values(cls, values: pd.Series) -> list[str]:
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


class PertAnnDataCatManager(CellxGeneAnnDataCatManager):
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

        self.PT_DEFAULT_VALUES = CellxGeneAnnDataCatManager.categoricals_defaults | {
            "cell_line": "unknown",
            "pert_target": "unknown",
        }

        self.PT_CATEGORICALS = CellxGeneAnnDataCatManager._get_cxg_categoricals() | {
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
    **kwargs,
):
    """Inspect values using a registry."""
    values = list(values)
    include_validated: list = []

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
            return {}
        return {"organism": organism or bt.settings.organism.name}
    return {}


def validate_categories(
    values: Iterable[str],
    field: FieldAttr,
    key: str,
    organism: str | None = None,
    source: Record | None = None,
    hint_print: str | None = None,
    curator: CatManager | None = None,
) -> tuple[bool, list[str]]:
    """Validate ontology terms using LaminDB registries.

    Args:
        values: The values to validate.
        field: The field attribute.
        key: The key referencing the slot in the DataFrame.
        organism: The organism name.
        source: The source record.
        standardize: Whether to standardize the values.
        hint_print: The hint to print that suggests fixing non-validated values.
    """
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
                warning_message += "\n    for remaining terms:\n"
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
    curator: CatManager | None = None,
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
            curator=curator,
            **kwargs,
        )
        validated &= is_val
        if len(non_val) > 0:
            non_validated[key] = non_val
    return validated, non_validated


def save_artifact(
    data: pd.DataFrame | ScverseDataStructures,
    *,
    fields: dict[str, FieldAttr] | dict[str, dict[str, FieldAttr]],
    index_field: FieldAttr | dict[str, FieldAttr] | None = None,
    description: str | None = None,
    organism: str | None = None,
    key: str | None = None,
    artifact: Artifact | None = None,
    revises: Artifact | None = None,
    run: Run | None = None,
    schema: Schema | None = None,
    **kwargs,
) -> Artifact:
    """Save all metadata with an Artifact.

    Args:
        data: The object to save.
        fields: A dictionary mapping obs_column to registry_field.
        index_field: The registry field to validate variables index against.
        description: A description of the artifact.
        organism: The organism name.
        key: A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a version family.
        artifact: A already registered artifact. Passing this will not save a new artifact from data.
        revises: Previous version of the artifact. Triggers a revision.
        run: The run that creates the artifact.
        schema: The Schema to associate with the Artifact.

    Returns:
        The saved Artifact.
    """
    from ..models.artifact import add_labels

    if artifact is None:
        if isinstance(data, pd.DataFrame):
            artifact = Artifact.from_df(
                data, description=description, key=key, revises=revises, run=run
            )
        elif isinstance(data, AnnData):
            artifact = Artifact.from_anndata(
                data, description=description, key=key, revises=revises, run=run
            )
        elif data_is_mudata(data):
            artifact = Artifact.from_mudata(
                data, description=description, key=key, revises=revises, run=run
            )
        elif data_is_spatialdata(data):
            artifact = Artifact.from_spatialdata(
                data, description=description, key=key, revises=revises, run=run
            )
        else:
            raise InvalidArgument(  # pragma: no cover
                "data must be one of pd.Dataframe, AnnData, MuData, SpatialData."
            )
    artifact.save()

    if organism is not None and index_field is not None:
        feature_kwargs = check_registry_organism(
            (
                list(index_field.values())[0].field.model
                if isinstance(index_field, dict)
                else index_field.field.model
            ),
            organism,
        )
    else:
        feature_kwargs = {}

    def _add_labels(
        data: pd.DataFrame | ScverseDataStructures,
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
            if not df[key].isna().all() and df[key].str.contains("|").any():
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

    match artifact.otype:
        case "DataFrame":
            artifact.features._add_set_from_df(field=index_field, **feature_kwargs)  # type: ignore
            _add_labels(
                data, artifact, fields, feature_ref_is_name=_ref_is_name(index_field)
            )
        case "AnnData":
            artifact.features._add_set_from_anndata(  # type: ignore
                var_field=index_field, **feature_kwargs
            )
            _add_labels(
                data, artifact, fields, feature_ref_is_name=_ref_is_name(index_field)
            )
        case "MuData":
            artifact.features._add_set_from_mudata(  # type: ignore
                var_fields=index_field, **feature_kwargs
            )
            for modality, modality_fields in fields.items():
                column_field_modality = index_field.get(modality)
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
        case "SpatialData":
            artifact.features._add_set_from_spatialdata(  # type: ignore
                sample_metadata_key=kwargs.get("sample_metadata_key", "sample"),
                var_fields=index_field,
                **feature_kwargs,
            )
            sample_metadata_key = kwargs.get("sample_metadata_key", "sample")
            for accessor, accessor_fields in fields.items():
                column_field = index_field.get(accessor)
                if accessor == sample_metadata_key:
                    _add_labels(
                        data.get_attrs(
                            key=sample_metadata_key, return_as="df", flatten=True
                        ),
                        artifact,
                        accessor_fields,
                        feature_ref_is_name=(
                            None if column_field is None else _ref_is_name(column_field)
                        ),
                    )
                else:
                    _add_labels(
                        data.tables[accessor],
                        artifact,
                        accessor_fields,
                        feature_ref_is_name=(
                            None if column_field is None else _ref_is_name(column_field)
                        ),
                    )
        case _:
            raise NotImplementedError  # pragma: no cover

    artifact.schema = schema
    artifact.save()

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
        kwargs: Additional keyword arguments to pass to the registry model to create new records.
    """
    from lamindb.models.save import save as ln_save

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
            save_ulabels_type(values, field=field, key=key)

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
    from ..models._from_values import _format_values

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


def save_ulabels_type(values: list[str], field: FieldAttr, key: str) -> None:
    """Save a parent label for the given labels."""
    registry = field.field.model
    assert registry == ULabel  # noqa: S101
    all_records = registry.filter(**{field.field.name: list(values)}).all()
    # so `tissue_type` becomes `TissueType`
    type_name = "".join([i.capitalize() for i in key.lower().split("_")])
    ulabel_type = registry.filter(name=type_name, is_type=True).one_or_none()
    if ulabel_type is None:
        ulabel_type = registry(name=type_name, is_type=True).save()
        logger.important(f"Created a ULabel type: {ulabel_type}")
    all_records.update(type=ulabel_type)


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


def _ref_is_name(field: FieldAttr | None) -> bool | None:
    """Check if the reference field is a name field."""
    from ..models.can_curate import get_name_field

    if field is not None:
        name_field = get_name_field(field.field.model)
        return field.field.name == name_field
    return None


# backward compat constructors ------------------


@classmethod  # type: ignore
def from_df(
    cls,
    df: pd.DataFrame,
    categoricals: dict[str, FieldAttr] | None = None,
    columns: FieldAttr = Feature.name,
    verbosity: str = "hint",
    organism: str | None = None,
) -> DataFrameCatManager:
    return DataFrameCatManager(
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
) -> AnnDataCatManager:
    return AnnDataCatManager(
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
    mdata: MuData | UPathStr,
    var_index: dict[str, dict[str, FieldAttr]],
    categoricals: dict[str, FieldAttr] | None = None,
    verbosity: str = "hint",
    organism: str | None = None,
) -> MuDataCatManager:
    return MuDataCatManager(
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
) -> TiledbsomaCatManager:
    return TiledbsomaCatManager(
        experiment_uri=experiment_uri,
        var_index=var_index,
        categoricals=categoricals,
        obs_columns=obs_columns,
        organism=organism,
        sources=sources,
    )


@classmethod  # type: ignore
def from_spatialdata(
    cls,
    sdata: SpatialData | UPathStr,
    var_index: dict[str, FieldAttr],
    categoricals: dict[str, dict[str, FieldAttr]] | None = None,
    organism: str | None = None,
    sources: dict[str, dict[str, Record]] | None = None,
    verbosity: str = "hint",
    *,
    sample_metadata_key: str = "sample",
):
    try:
        import spatialdata
    except ImportError as e:
        raise ImportError("Please install spatialdata: pip install spatialdata") from e

    return SpatialDataCatManager(
        sdata=sdata,
        var_index=var_index,
        categoricals=categoricals,
        verbosity=verbosity,
        organism=organism,
        sources=sources,
        sample_metadata_key=sample_metadata_key,
    )


CatManager.from_df = from_df  # type: ignore
CatManager.from_anndata = from_anndata  # type: ignore
CatManager.from_mudata = from_mudata  # type: ignore
CatManager.from_spatialdata = from_spatialdata  # type: ignore
CatManager.from_tiledbsoma = from_tiledbsoma  # type: ignore
