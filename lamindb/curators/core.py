"""Curator utilities.

.. autosummary::
   :toctree: .

   Curator
   SlotsCurator
   CatVector
   CatLookup
   DataFrameCatManager

"""

from __future__ import annotations

import copy
import re
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Callable

import lamindb_setup as ln_setup
import numpy as np
import pandas as pd
import pandera.pandas as pa
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args

from lamindb.base.types import FieldAttr  # noqa
from lamindb.models import (
    Artifact,
    Feature,
    Run,
    Schema,
    SQLRecord,
)
from lamindb.models._from_values import _format_values
from lamindb.models.artifact import (
    data_is_scversedatastructure,
    data_is_soma_experiment,
)
from lamindb.models.feature import parse_cat_dtype, parse_dtype

from ..errors import InvalidArgument, ValidationError

if TYPE_CHECKING:
    from typing import Any

    from anndata import AnnData
    from mudata import MuData
    from spatialdata import SpatialData
    from tiledbsoma._experiment import Experiment as SOMAExperiment

    from lamindb.core.types import ScverseDataStructures
    from lamindb.models.query_set import SQLRecordList


def strip_ansi_codes(text):
    # This pattern matches ANSI escape sequences
    ansi_pattern = re.compile(r"\x1b\[[0-9;]*m")
    return ansi_pattern.sub("", text)


class CatLookup:
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
        categoricals: list[Feature] | dict[str, FieldAttr],
        slots: dict[str, FieldAttr] = None,
        public: bool = False,
        sources: dict[str, SQLRecord] | None = None,
    ) -> None:
        slots = slots or {}
        if isinstance(categoricals, list):
            categoricals = {
                feature.name: parse_dtype(feature.dtype)[0]["field"]
                for feature in categoricals
            }
        self._categoricals = {**categoricals, **slots}
        self._public = public
        self._sources = sources

    def __getattr__(self, name):
        if name in self._categoricals:
            registry = self._categoricals[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public(source=self._sources.get(name)).lookup()
            else:
                return registry.lookup()
        raise AttributeError(
            f'"{self.__class__.__name__}" object has no attribute "{name}"'
        )

    def __getitem__(self, name):
        if name in self._categoricals:
            registry = self._categoricals[name].field.model
            if self._public and hasattr(registry, "public"):
                return registry.public(source=self._sources.get(name)).lookup()
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


SLOTS_DOCSTRING = """Access sub curators by slot."""


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
    """Curator base class.

    A `Curator` object makes it easy to validate, standardize & annotate datasets.

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
                self._dataset = self._dataset.load(is_run_input=False)
        self._schema: Schema | None = schema
        self._is_validated: bool = False

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
        pass  # pragma: no cover

    def __repr__(self) -> str:
        from lamin_utils import colors

        if self._schema is not None:
            # Schema might have different attributes
            if hasattr(self._schema, "name") and self._schema.name:
                schema_str = colors.italic(self._schema.name)
            elif hasattr(self._schema, "uid"):
                schema_str = colors.italic(f"uid={self._schema.uid}")
            elif hasattr(self._schema, "id"):
                schema_str = colors.italic(f"id={self._schema.id}")
            else:
                schema_str = colors.italic("unnamed")

            # Add schema type info if available
            if hasattr(self._schema, "otype") and self._schema.otype:
                schema_str += f" ({self._schema.otype})"
        else:
            schema_str = colors.warning("None")

        status_str = ""
        if self._is_validated:
            status_str = f", {colors.green('validated')}"
        else:
            status_str = f", {colors.yellow('unvalidated')}"

        cls_name = colors.green(self.__class__.__name__)

        # Get additional info based on curator type
        extra_info = ""
        if hasattr(self, "_slots") and self._slots:
            # For SlotsCurator and its subclasses
            slots_count = len(self._slots)
            if slots_count > 0:
                slot_names = list(self._slots.keys())
                if len(slot_names) <= 3:
                    extra_info = f", slots: {slot_names}"
                else:
                    extra_info = f", slots: [{', '.join(slot_names[:3])}... +{len(slot_names) - 3} more]"
        elif (
            cls_name == "DataFrameCurator"
            and hasattr(self, "cat")
            and hasattr(self.cat, "_categoricals")
        ):
            # For DataFrameCurator
            cat_count = len(getattr(self.cat, "_categoricals", []))
            if cat_count > 0:
                extra_info = f", categorical_features={cat_count}"

        artifact_info = ""
        if self._artifact is not None:
            artifact_uid = getattr(self._artifact, "uid", str(self._artifact))
            short_uid = (
                str(artifact_uid)[:8] + "..."
                if len(str(artifact_uid)) > 8
                else str(artifact_uid)
            )
            artifact_info = f", artifact: {colors.italic(short_uid)}"

        return (
            f"{cls_name}{artifact_info}(Schema: {schema_str}{extra_info}{status_str})"
        )


class SlotsCurator(Curator):
    """Curator for a dataset with slots.

    Args:
        dataset: The dataset to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    """

    def __init__(
        self,
        dataset: Artifact | ScverseDataStructures | SOMAExperiment,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        self._slots: dict[str, DataFrameCurator] = {}

        # used for multimodal data structures (not AnnData)
        # in form of {table/modality_key: var_field}
        self._var_fields: dict[str, FieldAttr] = {}
        # in form of {table/modality_key: categoricals}
        self._cat_vectors: dict[str, dict[str, CatVector]] = {}

    @property
    @doc_args(SLOTS_DOCSTRING)
    def slots(self) -> dict[str, DataFrameCurator]:
        """{}"""  # noqa: D415
        return self._slots

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}"""  # noqa: D415
        for slot, curator in self._slots.items():
            logger.info(f"validating slot {slot} ...")
            curator.validate()
        # set _is_validated to True as no slot raised an error
        self._is_validated = True

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

        if self._artifact is None:
            type_mapping = [
                (
                    lambda data: data_is_scversedatastructure(data, "AnnData"),
                    Artifact.from_anndata,
                ),
                (
                    lambda data: data_is_scversedatastructure(data, "MuData"),
                    Artifact.from_mudata,
                ),
                (
                    lambda data: data_is_scversedatastructure(data, "SpatialData"),
                    Artifact.from_spatialdata,
                ),
                (data_is_soma_experiment, Artifact.from_tiledbsoma),
            ]

            for type_check, factory in type_mapping:
                if type_check(self._dataset):
                    self._artifact = factory(  # type: ignore
                        self._dataset,
                        key=key,
                        description=description,
                        revises=revises,
                        run=run,
                    )
                    break

            self._artifact.schema = self._schema
            self._artifact.save()
        cat_vectors = {}
        for curator in self._slots.values():
            for key, cat_vector in curator.cat._cat_vectors.items():
                cat_vectors[key] = cat_vector
        return annotate_artifact(  # type: ignore
            self._artifact,
            curator=self,
            cat_vectors=cat_vectors,
        )


def is_list_of_type(value, expected_type):
    """Helper function to check if a value is either of expected_type or a list of that type, or a mix of both in a nested structure."""
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        # handle nested lists recursively
        return all(is_list_of_type(item, expected_type) for item in value)
    return isinstance(value, expected_type)


def check_dtype(expected_type) -> Callable:
    """Creates a check function for Pandera that validates a column's dtype.

    Supports both standard dtype checking and mixed list/single values for
    the same type. For example, a column with expected_type 'float' would
    also accept a mix of float values and lists of floats.

    Args:
        expected_type: String identifier for the expected type ('int', 'float', 'num', 'str')

    Returns:
        A function that checks if a series has the expected dtype or contains mixed types
    """

    def check_function(series):
        # first check if the series is entirely of the expected dtype (fast path)
        if expected_type == "int" and pd.api.types.is_integer_dtype(series.dtype):
            return True
        elif expected_type == "float" and pd.api.types.is_float_dtype(series.dtype):
            return True
        elif expected_type == "num" and pd.api.types.is_numeric_dtype(series.dtype):
            return True
        elif expected_type == "str" and pd.api.types.is_string_dtype(series.dtype):
            return True

        # if we're here, it might be a mixed column with object dtype
        # need to check each value individually
        if series.dtype == "object" and expected_type.startswith("list"):
            expected_type_member = expected_type.replace("list[", "").removesuffix("]")
            if expected_type_member == "int":
                return series.apply(lambda x: is_list_of_type(x, int)).all()
            elif expected_type_member == "float":
                return series.apply(lambda x: is_list_of_type(x, float)).all()
            elif expected_type_member == "num":
                # for numeric, accept either int or float
                return series.apply(lambda x: is_list_of_type(x, (int, float))).all()
            elif expected_type_member == "str" or expected_type_member.startswith(
                "cat["
            ):
                return series.apply(lambda x: is_list_of_type(x, str)).all()

        # if we get here, the validation failed
        return False

    return check_function


# this is also currently used as DictCurator
class DataFrameCurator(Curator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for `DataFrame`.

    Args:
        dataset: The DataFrame-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.
        slot: Indicate the slot in a composite curator for a composite data structure.

    Example:

        For simple example using a flexible schema, see :meth:`~lamindb.Artifact.from_df`.

        Here is an example that enforces a minimal set of columns in the dataframe.

        .. literalinclude:: scripts/curate_dataframe_minimal_errors.py
            :language: python

        Under-the-hood, this used the following schema.

        .. literalinclude:: scripts/define_mini_immuno_schema_flexible.py
            :language: python

        Valid features & labels were defined as:

        .. literalinclude:: scripts/define_mini_immuno_features_labels.py
            :language: python
    """

    def __init__(
        self,
        dataset: pd.DataFrame | Artifact,
        schema: Schema,
        slot: str | None = None,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        categoricals = []
        features = []
        feature_ids: set[int] = set()
        if schema.flexible:
            features += Feature.filter(name__in=self._dataset.keys()).list()
            feature_ids = {feature.id for feature in features}
        if schema.n > 0:
            if schema._index_feature_uid is not None:
                schema_features = [
                    feature
                    for feature in schema.members.list()
                    if feature.uid != schema._index_feature_uid  # type: ignore
                ]
            else:
                schema_features = schema.members.list()  # type: ignore
            if feature_ids:
                features.extend(
                    feature
                    for feature in schema_features
                    if feature.id not in feature_ids  # type: ignore
                )
            else:
                features.extend(schema_features)
        else:
            assert schema.itype is not None  # noqa: S101
        pandera_columns = {}
        if features or schema._index_feature_uid is not None:
            # populate features
            if schema.minimal_set:
                optional_feature_uids = set(schema.optionals.get_uids())
            for feature in features:
                if schema.minimal_set:
                    required = feature.uid not in optional_feature_uids
                else:
                    required = False
                # series.dtype is "object" if the column has lists types, e.g. [["a", "b"], ["a"], ["b"]]
                if feature.dtype in {"int", "float", "num"} or feature.dtype.startswith(
                    "list"
                ):
                    if isinstance(self._dataset, pd.DataFrame):
                        dtype = (
                            self._dataset[feature.name].dtype
                            if feature.name in self._dataset.keys()
                            else None
                        )
                    else:
                        dtype = None
                    pandera_columns[feature.name] = pa.Column(
                        dtype=None,
                        checks=pa.Check(
                            check_dtype(feature.dtype),
                            element_wise=False,
                            error=f"Column '{feature.name}' failed dtype check for '{feature.dtype}': got {dtype}",
                        ),
                        nullable=feature.nullable,
                        coerce=feature.coerce_dtype,
                        required=required,
                    )
                else:
                    pandera_dtype = (
                        feature.dtype
                        if not feature.dtype.startswith("cat")
                        else "category"
                    )
                    pandera_columns[feature.name] = pa.Column(
                        pandera_dtype,
                        nullable=feature.nullable,
                        coerce=feature.coerce_dtype,
                        required=required,
                    )
                if feature.dtype.startswith("cat") or feature.dtype.startswith(
                    "list[cat["
                ):
                    # validate categoricals if the column is required or if the column is present
                    if required or feature.name in self._dataset.keys():
                        categoricals.append(feature)
            if schema._index_feature_uid is not None:
                # in almost no case, an index should have a pandas.CategoricalDtype in a DataFrame
                # so, we're typing it as `str` here
                index = pa.Index(
                    schema.index.dtype
                    if not schema.index.dtype.startswith("cat")
                    else str
                )
            else:
                index = None
            self._pandera_schema = pa.DataFrameSchema(
                pandera_columns,
                coerce=schema.coerce_dtype,
                strict=schema.maximal_set,
                ordered=schema.ordered_set,
                index=index,
            )
        self._cat_manager = DataFrameCatManager(
            self._dataset,
            columns_field=parse_cat_dtype(schema.itype, is_itype=True)["field"],
            columns_names=pandera_columns.keys(),
            categoricals=categoricals,
            index=schema.index,
            slot=slot,
            maximal_set=schema.maximal_set,
        )

    @property
    @doc_args(CAT_MANAGER_DOCSTRING)
    def cat(self) -> DataFrameCatManager:
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
        self.cat.validate()

        if self.cat._is_validated:
            self._is_validated = True
        else:
            self._is_validated = False
            raise ValidationError(self.cat._validate_category_error_messages)

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}"""  # noqa: D415
        if self._schema.n > 0:
            try:
                # first validate through pandera
                self._pandera_schema.validate(self._dataset)
                # then validate lamindb categoricals
                self._cat_manager_validate()
            except pa.errors.SchemaError as err:
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
        if self._artifact is None:
            self._artifact = Artifact.from_df(
                self._dataset,
                key=key,
                description=description,
                revises=revises,
                run=run,
                format=".csv" if key.endswith(".csv") else None,
            )
            self._artifact.schema = self._schema
            self._artifact.save()
        return annotate_artifact(  # type: ignore
            self._artifact,
            cat_vectors=self.cat._cat_vectors,
        )


class AnnDataCurator(SlotsCurator):
    """Curator for `AnnData`.

    Args:
        dataset: The AnnData-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Example:

        .. literalinclude:: scripts/curate_anndata_flexible.py
            :language: python
            :caption: curate_anndata_flexible.py

    See Also:
        :meth:`~lamindb.Artifact.from_anndata`.
    """

    def __init__(
        self,
        dataset: AnnData | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_scversedatastructure(self._dataset, "AnnData"):
            raise InvalidArgument("dataset must be AnnData-like.")
        if schema.otype != "AnnData":
            raise InvalidArgument("Schema otype must be 'AnnData'.")
        self._slots = {
            slot: DataFrameCurator(
                (
                    getattr(self._dataset, slot.strip(".T")).T
                    if slot == "var.T"
                    or (
                        # backward compat
                        slot == "var"
                        and schema.slots["var"].itype not in {None, "Feature"}
                    )
                    else getattr(self._dataset, slot)
                ),
                slot_schema,
                slot=slot,
            )
            for slot, slot_schema in schema.slots.items()
            if slot in {"obs", "var", "var.T", "uns"}
        }
        if "var" in self._slots and schema.slots["var"].itype not in {None, "Feature"}:
            logger.warning(
                "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
            )
            self._slots["var"].cat._cat_vectors["var_index"] = self._slots[
                "var"
            ].cat._cat_vectors.pop("columns")
            self._slots["var"].cat._cat_vectors["var_index"]._key = "var_index"


def _assign_var_fields_categoricals_multimodal(
    modality: str | None,
    slot_type: str,
    slot: str,
    slot_schema: Schema,
    var_fields: dict[str, FieldAttr],
    cat_vectors: dict[str, dict[str, CatVector]],
    slots: dict[str, DataFrameCurator],
) -> None:
    """Assigns var_fields and categoricals for multimodal data curators."""
    if modality is not None:
        # Makes sure that all tables are present
        var_fields[modality] = None
        cat_vectors[modality] = {}

    if slot_type == "var":
        var_field = parse_cat_dtype(slot_schema.itype, is_itype=True)["field"]
        if modality is None:
            # This should rarely/never be used since tables should have different var fields
            var_fields[slot] = var_field  # pragma: no cover
        else:
            # Note that this is NOT nested since the nested key is always "var"
            var_fields[modality] = var_field
    else:
        obs_fields = slots[slot].cat._cat_vectors
        if modality is None:
            cat_vectors[slot] = obs_fields
        else:
            # Note that this is NOT nested since the nested key is always "obs"
            cat_vectors[modality] = obs_fields


class MuDataCurator(SlotsCurator):
    """Curator for `MuData`.

    Args:
        dataset: The MuData-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Example:

        .. literalinclude:: scripts/curate_mudata.py
            :language: python
            :caption: curate_mudata.py

    See Also:
        :meth:`~lamindb.Artifact.from_mudata`.
    """

    def __init__(
        self,
        dataset: MuData | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_scversedatastructure(self._dataset, "MuData"):
            raise InvalidArgument("dataset must be MuData-like.")
        if schema.otype != "MuData":
            raise InvalidArgument("Schema otype must be 'MuData'.")

        for slot, slot_schema in schema.slots.items():
            if ":" in slot:
                modality, modality_slot = slot.split(":")
                schema_dataset = self._dataset.__getitem__(modality)
            else:
                modality, modality_slot = None, slot
                schema_dataset = self._dataset
            if modality_slot == "var" and schema.slots[slot].itype not in {
                None,
                "Feature",
            }:
                logger.warning(
                    "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                )
            self._slots[slot] = DataFrameCurator(
                (
                    getattr(schema_dataset, modality_slot.rstrip(".T")).T
                    if modality_slot == "var.T"
                    or (
                        # backward compat
                        modality_slot == "var"
                        and schema.slots[slot].itype not in {None, "Feature"}
                    )
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
                cat_vectors=self._cat_vectors,
                slots=self._slots,
            )
        self._columns_field = self._var_fields


class SpatialDataCurator(SlotsCurator):
    """Curator for `SpatialData`.

    Args:
        dataset: The SpatialData-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Example:

        .. literalinclude:: scripts/curate_mudata.py
            :language: python
            :caption: curate_mudata.py

    See Also:
        :meth:`~lamindb.Artifact.from_spatialdata`.
    """

    def __init__(
        self,
        dataset: SpatialData | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_scversedatastructure(self._dataset, "SpatialData"):
            raise InvalidArgument("dataset must be SpatialData-like.")
        if schema.otype != "SpatialData":
            raise InvalidArgument("Schema otype must be 'SpatialData'.")

        for slot, slot_schema in schema.slots.items():
            split_result = slot.split(":")
            if (len(split_result) == 2 and split_result[0] == "table") or (
                len(split_result) == 3 and split_result[0] == "tables"
            ):
                if len(split_result) == 2:
                    table_key, sub_slot = split_result
                    logger.warning(
                        f"please prefix slot {slot} with 'tables:' going forward"
                    )
                else:
                    table_key, sub_slot = split_result[1], split_result[2]
                slot_object = self._dataset.tables.__getitem__(table_key)
                if sub_slot == "var" and schema.slots[slot].itype not in {
                    None,
                    "Feature",
                }:
                    logger.warning(
                        "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                    )
                data_object = (
                    getattr(slot_object, sub_slot.rstrip(".T")).T
                    if sub_slot == "var.T"
                    or (
                        # backward compat
                        sub_slot == "var"
                        and schema.slots[slot].itype not in {None, "Feature"}
                    )
                    else getattr(slot_object, sub_slot)
                )
            elif len(split_result) == 1 or (
                len(split_result) > 1 and split_result[0] == "attrs"
            ):
                table_key = None
                if len(split_result) == 1:
                    if split_result[0] != "attrs":
                        logger.warning(
                            f"please prefix slot {slot} with 'attrs:' going forward"
                        )
                        sub_slot = slot
                        data_object = self._dataset.attrs[slot]
                    else:
                        sub_slot = "attrs"
                        data_object = self._dataset.attrs
                elif len(split_result) == 2:
                    sub_slot = split_result[1]
                    data_object = self._dataset.attrs[split_result[1]]
                data_object = pd.DataFrame([data_object])
            self._slots[slot] = DataFrameCurator(data_object, slot_schema, slot)
            _assign_var_fields_categoricals_multimodal(
                modality=table_key,
                slot_type=sub_slot,
                slot=slot,
                slot_schema=slot_schema,
                var_fields=self._var_fields,
                cat_vectors=self._cat_vectors,
                slots=self._slots,
            )
        self._columns_field = self._var_fields


class TiledbsomaExperimentCurator(SlotsCurator):
    """Curator for `TileDB-SOMA`.

    Args:
        dataset: The `tiledbsoma.Experiment` object.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Example:

        .. literalinclude:: scripts/curate_soma_experiment.py
            :language: python
            :caption: curate_soma_experiment.py

    See Also:
        :meth:`~lamindb.Artifact.from_tiledbsoma`.
    """

    def __init__(
        self,
        dataset: SOMAExperiment | Artifact,
        schema: Schema,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema)
        if not data_is_soma_experiment(self._dataset):
            raise InvalidArgument("dataset must be SOMAExperiment-like.")
        if schema.otype != "tiledbsoma":
            raise InvalidArgument("Schema otype must be 'tiledbsoma'.")

        for slot, slot_schema in schema.slots.items():
            if slot.startswith("ms:"):
                ms, modality_slot = slot.split(":")
                schema_dataset = (
                    self._dataset.ms[modality_slot.removesuffix(".T")]
                    .var.read()
                    .concat()
                    .to_pandas()
                    .drop("soma_joinid", axis=1, errors="ignore")
                )

                self._slots[slot] = DataFrameCurator(
                    (
                        schema_dataset.T
                        if modality_slot == "var.T"
                        or (
                            # backward compat
                            modality_slot == "var"
                            and schema.slots[slot].itype not in {None, "Feature"}
                        )
                        else schema_dataset
                    ),
                    slot_schema,
                )
            else:
                # global Experiment obs slot
                _ms, modality_slot = None, slot
                schema_dataset = (
                    self._dataset.obs.read()
                    .concat()
                    .to_pandas()
                    .drop(["soma_joinid", "obs_id"], axis=1, errors="ignore")
                )
                self._slots[slot] = DataFrameCurator(
                    schema_dataset,
                    slot_schema,
                )

            if modality_slot == "var" and schema.slots[slot].itype not in {
                None,
                "Feature",
            }:
                logger.warning(
                    "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                )

            _assign_var_fields_categoricals_multimodal(
                modality=slot,  # not using "ms" here as it would always be the same for all modalities
                slot_type=modality_slot,
                slot=slot,
                slot_schema=slot_schema,
                var_fields=self._var_fields,
                cat_vectors=self._cat_vectors,
                slots=self._slots,
            )
        self._columns_field = self._var_fields


class CatVector:
    """Vector with categorical values."""

    def __init__(
        self,
        values_getter: Callable
        | Iterable[str],  # A callable or iterable that returns the values to validate.
        field: FieldAttr,  # The field to validate against.
        key: str,  # The name of the vector to validate. Only used for logging.
        values_setter: Callable | None = None,  # A callable that sets the values.
        source: SQLRecord | None = None,  # The ontology source to validate against.
        feature: Feature | None = None,
        cat_manager: DataFrameCatManager | None = None,
        subtype_str: str = "",
        maximal_set: bool = True,  # whether unvalidated categoricals cause validation failure.
    ) -> None:
        self._values_getter = values_getter
        self._values_setter = values_setter
        self._field = field
        self._key = key
        self._source = source
        self._organism = None
        self._validated: None | list[str] = None
        self._non_validated: None | list[str] = None
        self._synonyms: None | dict[str, str] = None
        self._subtype_str = subtype_str
        self._subtype_query_set = None
        self._cat_manager = cat_manager
        self.feature = feature
        self.records = None
        self._maximal_set = maximal_set
        if hasattr(field.field.model, "_name_field"):
            label_ref_is_name = field.field.name == field.field.model._name_field
        else:
            label_ref_is_name = field.field.name == "name"
        self.label_ref_is_name = label_ref_is_name

    @property
    def values(self):
        """Get the current values using the getter function."""
        if callable(self._values_getter):
            return self._values_getter()
        return self._values_getter

    @values.setter
    def values(self, new_values):
        """Set new values using the setter function if available."""
        if callable(self._values_setter):
            self._values_setter(new_values)
        else:
            # If values_getter is not callable, it's a direct reference we can update
            self._values_getter = new_values

    @property
    def is_validated(self) -> bool:
        """Whether the vector is validated."""
        # if nothing was validated, something likely is fundamentally wrong
        # should probably add a setting `at_least_one_validated`
        result = True
        if len(self.values) > 0 and len(self.values) == len(self._non_validated):
            result = False
        # len(self._non_validated) != 0
        #     if maximal_set is True, return False
        #     if maximal_set is False, return True
        # len(self._non_validated) == 0
        #     return True
        if len(self._non_validated) != 0:
            if self._maximal_set:
                result = False
        return result

    def _replace_synonyms(self) -> list[str]:
        """Replace synonyms in the vector with standardized values."""

        def process_value(value, syn_mapper):
            """Helper function to process values recursively."""
            if isinstance(value, list):
                # Handle list - recursively process each item
                return [process_value(item, syn_mapper) for item in value]
            else:
                # Handle single value
                return syn_mapper.get(value, value)

        syn_mapper = self._synonyms
        # replace the values in df
        std_values = self.values.map(
            lambda unstd_val: process_value(unstd_val, syn_mapper)
        )
        # remove the standardized values from self.non_validated
        non_validated = [i for i in self._non_validated if i not in syn_mapper]
        if len(non_validated) == 0:
            self._non_validated = []
        else:
            self._non_validated = non_validated  # type: ignore
        # logging
        n = len(syn_mapper)
        if n > 0:
            syn_mapper_print = _format_values(
                [f'"{k}" → "{v}"' for k, v in syn_mapper.items()], sep=""
            )
            s = "s" if n > 1 else ""
            logger.success(
                f'standardized {n} synonym{s} in "{self._key}": {colors.green(syn_mapper_print)}'
            )
        return std_values

    def __repr__(self) -> str:
        if self._non_validated is None:
            status = "unvalidated"
        else:
            status = (
                "validated"
                if len(self._non_validated) == 0
                else f"non-validated ({len(self._non_validated)})"
            )

        field_name = getattr(self._field, "name", str(self._field))
        values_count = len(self.values) if hasattr(self.values, "__len__") else "?"
        return f"CatVector(key='{self._key}', field='{field_name}', values={values_count}, {status})"

    def _add_validated(self) -> tuple[list, list]:
        """Save features or labels records in the default instance."""
        from lamindb.models.save import save as ln_save

        registry = self._field.field.model
        field_name = self._field.field.name
        model_field = registry.__get_name_with_module__()
        filter_kwargs = get_current_filter_kwargs(
            registry, {"organism": self._organism, "source": self._source}
        )
        values = [
            i
            for i in self.values
            if (isinstance(i, str) and i)
            or (isinstance(i, list) and i)
            or (isinstance(i, np.ndarray) and i.size > 0)
        ]
        if not values:
            return [], []

        # if a value is a list, we need to flatten it
        str_values = _flatten_unique(values)

        # inspect the default instance and save validated records from public
        if (
            self._subtype_str != "" and "__" not in self._subtype_str
        ):  # not for general filter expressions
            related_name = registry._meta.get_field("type").remote_field.related_name
            self._subtype_query_set = getattr(
                registry.get(name=self._subtype_str), related_name
            ).all()
            values_array = np.array(str_values)
            validated_mask = self._subtype_query_set.validate(  # type: ignore
                values_array, field=self._field, **filter_kwargs, mute=True
            )
            validated_labels, non_validated_labels = (
                values_array[validated_mask],
                values_array[~validated_mask],
            )
            records = registry.from_values(
                validated_labels, field=self._field, **filter_kwargs, mute=True
            )
        else:
            existing_and_public_records = registry.from_values(
                str_values, field=self._field, **filter_kwargs, mute=True
            )
            existing_and_public_labels = [
                getattr(r, field_name) for r in existing_and_public_records
            ]
            # public records that are not already in the database
            public_records = [r for r in existing_and_public_records if r._state.adding]
            # here we check to only save the public records if they are from the specified source
            # we check the uid because r.source and source can be from different instances
            if self._source:
                public_records = [
                    r for r in public_records if r.source.uid == self._source.uid
                ]
            if len(public_records) > 0:
                logger.info(f"saving validated records of '{self._key}'")
                ln_save(public_records)
                labels_saved_public = [getattr(r, field_name) for r in public_records]
                # log the saved public labels
                # the term "transferred" stresses that this is always in the context of transferring
                # labels from a public ontology or a different instance to the present instance
                if len(labels_saved_public) > 0:
                    s = "s" if len(labels_saved_public) > 1 else ""
                    logger.success(
                        f'added {len(labels_saved_public)} record{s} {colors.green("from_public")} with {model_field} for "{self._key}": {_format_values(labels_saved_public)}'
                    )
                    # non-validated records from the default instance
            non_validated_labels = [
                i for i in str_values if i not in existing_and_public_labels
            ]
            validated_labels = existing_and_public_labels
            records = existing_and_public_records

        self.records = records
        # validated, non-validated
        return validated_labels, non_validated_labels

    def _add_new(
        self,
        values: list[str],
        df: pd.DataFrame | None = None,  # remove when all users use schema
        dtype: str | None = None,
        **create_kwargs,
    ) -> None:
        """Add new labels to the registry."""
        from lamindb.models.save import save as ln_save

        registry = self._field.field.model
        field_name = self._field.field.name
        non_validated_records: SQLRecordList[Any] = []  # type: ignore
        if df is not None and registry == Feature:
            nonval_columns = Feature.inspect(df.columns, mute=True).non_validated
            non_validated_records = Feature.from_df(df.loc[:, nonval_columns])
        else:
            if (
                self._organism
                and hasattr(registry, "organism")
                and registry._meta.get_field("organism").is_relation
            ):
                # make sure organism record is saved to the current instance
                create_kwargs["organism"] = _save_organism(name=self._organism)

            for value in values:
                init_kwargs = {field_name: value}
                if registry == Feature:
                    init_kwargs["dtype"] = "cat" if dtype is None else dtype
                non_validated_records.append(registry(**init_kwargs, **create_kwargs))
        if len(non_validated_records) > 0:
            ln_save(non_validated_records)
            model_field = colors.italic(registry.__get_name_with_module__())
            s = "s" if len(values) > 1 else ""
            logger.success(
                f'added {len(values)} record{s} with {model_field} for "{self._key}": {_format_values(values)}'
            )

    def _validate(
        self,
        values: list[str],
    ) -> tuple[list[str], dict]:
        """Validate ontology terms using LaminDB registries."""
        registry = self._field.field.model
        field_name = self._field.field.name
        model_field = f"{registry.__name__}.{field_name}"

        kwargs_current = get_current_filter_kwargs(
            registry, {"organism": self._organism, "source": self._source}
        )

        # inspect values from the default instance, excluding public
        registry_or_queryset = registry
        if self._subtype_query_set is not None:
            registry_or_queryset = self._subtype_query_set
        inspect_result = registry_or_queryset.inspect(
            values, field=self._field, mute=True, from_source=False, **kwargs_current
        )
        non_validated = inspect_result.non_validated
        syn_mapper = inspect_result.synonyms_mapper

        # inspect the non-validated values from public (BioRecord only)
        values_validated = []
        if hasattr(registry, "public"):
            public_records = registry.from_values(
                non_validated,
                field=self._field,
                mute=True,
                **kwargs_current,
            )
            values_validated += [getattr(r, field_name) for r in public_records]

        # logging messages
        if self._cat_manager is not None:
            slot = self._cat_manager._slot
        else:
            slot = None
        in_slot = f" in slot '{slot}'" if slot is not None else ""
        slot_prefix = f".slots['{slot}']" if slot is not None else ""
        non_validated_hint_print = (
            f"curator{slot_prefix}.cat.add_new_from('{self._key}')"
        )
        non_validated = [i for i in non_validated if i not in values_validated]
        n_non_validated = len(non_validated)
        if n_non_validated == 0:
            logger.success(
                f'"{self._key}" is validated against {colors.italic(model_field)}'
            )
            return [], {}
        else:
            s = "" if n_non_validated == 1 else "s"
            print_values = _format_values(non_validated)
            warning_message = f"{colors.red(f'{n_non_validated} term{s}')} not validated in feature '{self._key}'{in_slot}: {colors.red(print_values)}\n"
            if syn_mapper:
                s = "" if len(syn_mapper) == 1 else "s"
                syn_mapper_print = _format_values(
                    [f'"{k}" → "{v}"' for k, v in syn_mapper.items()], sep=""
                )
                hint_msg = f'.standardize("{self._key}")'
                warning_message += f"    {colors.yellow(f'{len(syn_mapper)} synonym{s}')} found: {colors.yellow(syn_mapper_print)}\n    → curate synonyms via: {colors.cyan(hint_msg)}"
            if n_non_validated > len(syn_mapper):
                if syn_mapper:
                    warning_message += "\n    for remaining terms:\n"
                warning_message += f"    → fix typos, remove non-existent values, or save terms via: {colors.cyan(non_validated_hint_print)}"
                if self._subtype_query_set is not None:
                    warning_message += f"\n    → a valid label for subtype '{self._subtype_str}' has to be one of {self._subtype_query_set.list('name')}"
            logger.info(f'mapping "{self._key}" on {colors.italic(model_field)}')
            logger.warning(warning_message)
            if self._cat_manager is not None:
                self._cat_manager._validate_category_error_messages = strip_ansi_codes(
                    warning_message
                )
            return non_validated, syn_mapper

    def validate(self) -> None:
        """Validate the vector."""
        # add source-validated values to the registry
        self._validated, self._non_validated = self._add_validated()
        self._non_validated, self._synonyms = self._validate(values=self._non_validated)

        # always register new Features if they are columns
        if self._key == "columns" and self._field == Feature.name:
            self.add_new()

    def standardize(self) -> None:
        """Standardize the vector."""
        registry = self._field.field.model
        if not hasattr(registry, "standardize"):
            return self.values
        if self._synonyms is None:
            self.validate()
        # get standardized values
        std_values = self._replace_synonyms()
        # update non_validated values
        self._non_validated = [
            i for i in self._non_validated if i not in self._synonyms.keys()
        ]
        # remove synonyms since they are now standardized
        self._synonyms = {}
        # update the values with the standardized values
        self.values = std_values

    def add_new(self, **create_kwargs) -> None:
        """Add new values to the registry."""
        if self._non_validated is None:
            self.validate()
        if len(self._synonyms) > 0:
            # raise error because .standardize modifies the input dataset
            raise ValidationError(
                "Please run `.standardize()` before adding new values."
            )
        self._add_new(
            values=self._non_validated,
            **create_kwargs,
        )
        # remove the non_validated values since they are now registered
        self._non_validated = []


class DataFrameCatManager:
    """Manage categoricals by updating registries.

    This class is accessible from within a `DataFrameCurator` via the `.cat` attribute.

    If you find non-validated values, you have two options:

    - new values found in the data can be registered via `DataFrameCurator.cat.add_new_from()` :meth:`~lamindb.curators.core.DataFrameCatManager.add_new_from`
    - non-validated values can be accessed via `DataFrameCurator.cat.add_new_from()` :meth:`~lamindb.curators.core.DataFrameCatManager.non_validated` and addressed manually
    """

    def __init__(
        self,
        df: pd.DataFrame | Artifact,
        columns_field: FieldAttr = Feature.name,
        columns_names: Iterable[str] | None = None,
        categoricals: list[Feature] | None = None,
        sources: dict[str, SQLRecord] | None = None,
        index: Feature | None = None,
        slot: str | None = None,
        maximal_set: bool = False,
    ) -> None:
        self._non_validated = None
        self._index = index
        self._artifact: Artifact = None  # pass the dataset as an artifact
        self._dataset: Any = df  # pass the dataset as a UPathStr or data object
        if isinstance(self._dataset, Artifact):
            self._artifact = self._dataset
            self._dataset = self._dataset.load(is_run_input=False)
        self._is_validated: bool = False
        self._categoricals = categoricals or []
        self._non_validated = None
        self._sources = sources or {}
        self._columns_field = columns_field
        self._validate_category_error_messages: str = ""
        self._cat_vectors: dict[str, CatVector] = {}
        self._slot = slot
        self._maximal_set = maximal_set

        if columns_names is None:
            columns_names = []
        if columns_field == Feature.name:
            self._cat_vectors["columns"] = CatVector(
                values_getter=columns_names,
                field=columns_field,
                key="columns" if isinstance(self._dataset, pd.DataFrame) else "keys",
                source=self._sources.get("columns"),
                cat_manager=self,
                maximal_set=self._maximal_set,
            )
        else:
            self._cat_vectors["columns"] = CatVector(
                values_getter=lambda: self._dataset.columns,  # lambda ensures the inplace update
                values_setter=lambda new_values: setattr(
                    self._dataset, "columns", pd.Index(new_values)
                ),
                field=columns_field,
                key="columns",
                source=self._sources.get("columns"),
                cat_manager=self,
                maximal_set=self._maximal_set,
            )
        for feature in self._categoricals:
            result = parse_dtype(feature.dtype)[
                0
            ]  # TODO: support composite dtypes for categoricals
            key = feature.name
            field = result["field"]
            subtype_str = result["subtype_str"]
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
                feature=feature,
                cat_manager=self,
                subtype_str=subtype_str,
            )
        if index is not None and index.dtype.startswith("cat"):
            result = parse_dtype(index.dtype)[0]
            field = result["field"]
            key = "index"
            self._cat_vectors[key] = CatVector(
                values_getter=self._dataset.index,
                field=field,
                key=key,
                feature=index,
                cat_manager=self,
            )

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
    def categoricals(self) -> list[Feature]:
        """The categorical features."""
        return self._categoricals

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
        for key, cat_vector in self._cat_vectors.items():
            logger.info(f"validating vector {key}")
            cat_vector.validate()
            validated &= cat_vector.is_validated
        self._is_validated = validated
        self._non_validated = {}  # type: ignore

        if self._index is not None:
            # cat_vector.validate() populates validated labels
            # the index should become part of the feature set corresponding to the dataframe
            if self._cat_vectors["columns"].records is not None:
                self._cat_vectors["columns"].records.insert(0, self._index)  # type: ignore
            else:
                self._cat_vectors["columns"].records = [self._index]  # type: ignore

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


def get_current_filter_kwargs(registry: type[SQLRecord], kwargs: dict) -> dict:
    """Make sure the source and organism are saved in the same database as the registry."""
    db = registry.filter().db
    source = kwargs.get("source")
    organism = kwargs.get("organism")
    filter_kwargs = kwargs.copy()

    if isinstance(organism, SQLRecord) and organism._state.db != "default":
        if db is None or db == "default":
            organism_default = copy.copy(organism)
            # save the organism record in the default database
            organism_default.save()
            filter_kwargs["organism"] = organism_default
    if isinstance(source, SQLRecord) and source._state.db != "default":
        if db is None or db == "default":
            source_default = copy.copy(source)
            # save the source record in the default database
            source_default.save()
            filter_kwargs["source"] = source_default

    return filter_kwargs


def get_organism_kwargs(
    field: FieldAttr, organism: str | None = None, values: Any = None
) -> dict[str, str]:
    """Check if a registry needs an organism and return the organism name."""
    registry = field.field.model
    if registry.__base__.__name__ == "BioRecord":
        import bionty as bt
        from bionty._organism import is_organism_required

        from ..models._from_values import get_organism_record_from_field

        if is_organism_required(registry):
            if organism is not None or bt.settings.organism is not None:
                return {"organism": organism or bt.settings.organism.name}
            else:
                organism_record = get_organism_record_from_field(
                    field, organism=organism, values=values
                )
                if organism_record is not None:
                    return {"organism": organism_record.name}
    return {}


def annotate_artifact(
    artifact: Artifact,
    *,
    curator: AnnDataCurator | SlotsCurator | None = None,
    cat_vectors: dict[str, CatVector] | None = None,
) -> Artifact:
    from .. import settings
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
        if len(cat_vector.records) > settings.annotation.n_max_records:
            logger.important(
                f"not annotating with {len(cat_vector.records)} labels for feature {key} as it exceeds {settings.annotation.n_max_records} (ln.settings.annotation.n_max_records)"
            )
            continue
        add_labels(
            artifact,
            records=cat_vector.records,
            feature=cat_vector.feature,
            feature_ref_is_name=None,  # do not need anymore
            label_ref_is_name=cat_vector.label_ref_is_name,
            from_curator=True,
        )

    # annotate with inferred schemas aka feature sets
    if artifact.otype == "DataFrame":
        features = cat_vectors["columns"].records
        if features is not None:
            feature_set = Schema(
                features=features, coerce_dtype=artifact.schema.coerce_dtype
            )  # TODO: add more defaults from validating schema
            if (
                feature_set._state.adding
                and len(features) > settings.annotation.n_max_records
            ):
                logger.important(
                    f"not annotating with {len(features)} features as it exceeds {settings.annotation.n_max_records} (ln.settings.annotation.n_max_records)"
                )
                itype = parse_cat_dtype(artifact.schema.itype, is_itype=True)["field"]
                feature_set = Schema(itype=itype, n=len(features))
            artifact.feature_sets.add(
                feature_set.save(), through_defaults={"slot": "columns"}
            )
    else:
        for slot, slot_curator in curator._slots.items():
            # var_index is backward compat (2025-05-01)
            name = (
                "var_index"
                if (slot == "var" and "var_index" in slot_curator.cat._cat_vectors)
                else "columns"
            )
            features = slot_curator.cat._cat_vectors[name].records
            if features is None:
                logger.warning(f"no features found for slot {slot}")
                continue
            itype = parse_cat_dtype(artifact.schema.slots[slot].itype, is_itype=True)[
                "field"
            ]
            feature_set = Schema(features=features, itype=itype)
            if (
                feature_set._state.adding
                and len(features) > settings.annotation.n_max_records
            ):
                logger.important(
                    f"not annotating with {len(features)} features for slot {slot} as it exceeds {settings.annotation.n_max_records} (ln.settings.annotation.n_max_records)"
                )
                feature_set = Schema(itype=itype, n=len(features))
            artifact.feature_sets.add(
                feature_set.save(), through_defaults={"slot": slot}
            )

    slug = ln_setup.settings.instance.slug
    if ln_setup.settings.instance.is_remote:  # pdagma: no cover
        logger.important(f"go to https://lamin.ai/{slug}/artifact/{artifact.uid}")
    return artifact


def _flatten_unique(series: pd.Series[list[Any] | Any]) -> list[Any]:
    """Flatten a Pandas series containing lists or single items into a unique list of elements.

    The order of elements in the result list preserves the order they first appear in the input series.
    """
    # Use dict.fromkeys to preserve order while ensuring uniqueness
    result: dict = {}

    for item in series:
        if isinstance(item, list | np.ndarray):
            # Add each element to the dict (only first occurrence is kept)
            for element in item:
                result[element] = None
        else:
            result[item] = None

    # Return the keys as a list, preserving order
    return list(result.keys())


def _save_organism(name: str):
    """Save an organism record."""
    import bionty as bt

    organism = bt.Organism.filter(name=name).one_or_none()
    if organism is None:
        organism = bt.Organism.from_source(name=name)
        if organism is None:
            raise ValidationError(
                f'Organism "{name}" not found from public reference\n'
                f'      → please save it from a different source: bt.Organism.from_source(name="{name}", source).save()'
                f'      → or manually save it without source: bt.Organism(name="{name}").save()'
            )
        organism.save()
    return organism
