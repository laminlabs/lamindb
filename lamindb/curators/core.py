"""Curator utilities.

.. autoclass:: Curator
.. autoclass:: SlotsCurator
.. autoclass:: ComponentCurator
.. autoclass:: CatVector
.. autoclass:: CatLookup
.. autoclass:: DataFrameCatManager

"""

from __future__ import annotations

import copy
import re
from typing import TYPE_CHECKING, Any, Callable

import lamindb_setup as ln_setup
import numpy as np
import pandas as pd
import pandera.pandas as pandera
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core.upath import LocalPathClasses

from lamindb.base.dtypes import check_dtype
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
from lamindb.models.feature import (
    parse_cat_dtype,
    parse_dtype,
    parse_filter_string,
    resolve_relation_filters,
)

from ..errors import InvalidArgument, ValidationError

if TYPE_CHECKING:
    from collections.abc import Iterable
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

SLOTS_DETAILS_DOCSTRING = """Uses **slots** to specify which component contains which schema. Slots are keys that identify where features are stored within composite data structures."""

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
        - :class:`~lamindb.curators.TiledbsomaExperimentCurator`
    """

    def __init__(
        self, dataset: Any, schema: Schema, *, features: dict[str, Any] | None = None
    ) -> None:
        if not isinstance(schema, Schema):
            raise InvalidArgument("schema argument must be a Schema record.")
        if schema.pk is None:
            raise ValueError(
                "Schema must be saved before curation. Please save it using '.save()'."
            )
        self._artifact: Artifact | None = None
        self._dataset: Any = None
        # self._dataset is set below, it is opened or loaded if dataset is an Artifact
        if isinstance(dataset, Artifact):
            self._artifact = dataset
            if self._artifact.otype in {
                "DataFrame",
                "AnnData",
                "MuData",
                "SpatialData",
            }:
                if (
                    not isinstance(self._artifact.path, LocalPathClasses)
                    and self._artifact.otype == "AnnData"
                ):
                    try:
                        self._dataset = self._artifact.open(mode="r")
                        logger.important(
                            "opened remote artifact for streaming during validation"
                        )
                    except Exception as e:
                        logger.warning(
                            f"unable to open remote AnnData Artifact: {e}, falling back to loading into memory"
                        )
                if self._dataset is None:
                    logger.important("loading artifact into memory for validation")
                    self._dataset = self._artifact.load(is_run_input=False)
            else:
                raise InvalidArgument(
                    f"Cannot load or open artifact of this type: {self._artifact}"
                )
        else:
            self._dataset = dataset
        self._schema: Schema = schema
        self._external_features: dict[str, Any] = features or {}
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
            artifact_info = f", artifact: {colors.italic(self._artifact.uid)}"

        return (
            f"{cls_name}{artifact_info}(Schema: {schema_str}{extra_info}{status_str})"
        )


@doc_args(SLOTS_DETAILS_DOCSTRING)
class SlotsCurator(Curator):
    """Curator for a dataset with slots.

    {}

    Args:
        dataset: The dataset to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.
    """

    def __init__(
        self,
        dataset: Artifact | ScverseDataStructures | SOMAExperiment,
        schema: Schema,
        *,
        features: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(dataset=dataset, schema=schema, features=features)
        self._slots: dict[str, ComponentCurator] = {}

        # used for multimodal data structures (not AnnData)
        # in form of {table/modality_key: var_field}
        self._var_fields: dict[str, FieldAttr] = {}
        # in form of {table/modality_key: categoricals}
        self._cat_vectors: dict[str, dict[str, CatVector]] = {}

    @property
    @doc_args(SLOTS_DOCSTRING)
    def slots(self) -> dict[str, ComponentCurator]:
        """{}"""  # noqa: D415
        return self._slots

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}"""  # noqa: D415
        if self._external_features is not None and "__external__" in self._schema.slots:
            validation_schema = self._schema.slots["__external__"]
            ExperimentalDictCurator(
                self._external_features, validation_schema
            ).validate()
        for slot, curator in self._slots.items():
            logger.debug(f"validating slot {slot} ...")
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
                    lambda dataset: isinstance(dataset, pd.DataFrame),
                    Artifact.from_dataframe,
                ),
                (
                    lambda dataset: data_is_scversedatastructure(dataset, "AnnData"),
                    Artifact.from_anndata,
                ),
                (
                    lambda dataset: data_is_scversedatastructure(dataset, "MuData"),
                    Artifact.from_mudata,
                ),
                (
                    lambda dataset: data_is_scversedatastructure(
                        dataset, "SpatialData"
                    ),
                    Artifact.from_spatialdata,
                ),
                (data_is_soma_experiment, Artifact.from_tiledbsoma),
            ]

            for type_check, af_constructor in type_mapping:
                if type_check(self._dataset):
                    # not all Artifact constructors support `features` yet
                    self._artifact = af_constructor(  # type: ignore
                        self._dataset,
                        key=key,
                        description=description,
                        revises=revises,
                        run=run,
                    )
                    self._artifact._external_features = self._external_features
                    break

        cat_vectors = {}
        for curator in self._slots.values():
            for key, cat_vector in curator.cat._cat_vectors.items():
                cat_vectors[key] = cat_vector

        self._artifact.schema = self._schema
        if self._external_features:
            assert self._artifact._external_features == self._external_features
        self._artifact.save()
        return annotate_artifact(  # type: ignore
            self._artifact,
            curator=self,
            cat_vectors=cat_vectors,
        )


def convert_dict_to_dataframe_for_validation(d: dict, schema: Schema) -> pd.DataFrame:
    """Convert a dictionary to a DataFrame for validation against a schema."""
    df = pd.DataFrame([d])
    for feature in schema.members:
        # we cannot cast a `list[cat[...]]]` to categorical because lists are not hashable
        if feature.dtype.startswith("cat"):
            if feature.name in df.columns:
                df[feature.name] = pd.Categorical(df[feature.name])
    return df


# For more context, read https://laminlabs.slack.com/archives/C07DB677JF6/p1753994077716099 and
# https://www.notion.so/laminlabs/Add-a-DictCurator-2422aeaa55e180b9a513f91d13970836
class ComponentCurator(Curator):
    """Curator for `DataFrame`.

    Provides all key functionality to validate Pandas DataFrames.
    This class is not user facing unlike :class:`~lamindb.curators.DataFrameCurator` which extends this
    class with functionality to validate the `attrs` slot.

    Args:
        dataset: The DataFrame-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.
        slot: Indicate the slot in a composite curator for a composite data structure.
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
            features += Feature.filter(name__in=self._dataset.keys()).to_list()
            feature_ids = {feature.id for feature in features}

        if schema.n > 0:
            if schema._index_feature_uid is not None:
                schema_features = [
                    feature
                    for feature in schema.members.to_list()
                    if feature.uid != schema._index_feature_uid  # type: ignore
                ]
            else:
                schema_features = schema.members.to_list()  # type: ignore
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
        self._pandera_schema = None
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
                if feature.dtype.startswith("list[cat"):
                    pandera_columns[feature.name] = pandera.Column(
                        dtype=None,
                        checks=pandera.Check(
                            check_dtype("list", feature.nullable),
                            element_wise=False,
                            error=f"Column '{feature.name}' failed dtype check for '{feature.dtype}'",
                        ),
                        nullable=feature.nullable,
                        coerce=feature.coerce_dtype,
                        required=required,
                    )
                elif feature.dtype in {
                    "int",
                    "float",
                    "num",
                    "path",
                } or feature.dtype.startswith("list"):
                    if isinstance(self._dataset, pd.DataFrame):
                        dtype = (
                            self._dataset[feature.name].dtype
                            if feature.name in self._dataset.keys()
                            else None
                        )
                    else:
                        dtype = None
                    pandera_columns[feature.name] = pandera.Column(
                        dtype=None,
                        checks=pandera.Check(
                            check_dtype(feature.dtype, feature.nullable),
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
                    pandera_columns[feature.name] = pandera.Column(
                        pandera_dtype,
                        nullable=feature.nullable,
                        coerce=feature.coerce_dtype,
                        required=required,
                    )
                if feature.dtype.startswith("cat") or feature.dtype.startswith(
                    "list[cat["
                ):
                    # validate categoricals if the column is required or if the column is present
                    # but exclude the index feature from column categoricals
                    if (required or feature.name in self._dataset.keys()) and (
                        schema._index_feature_uid is None
                        or feature.uid != schema._index_feature_uid
                    ):
                        categoricals.append(feature)
            # in almost no case, an index should have a pandas.CategoricalDtype in a DataFrame
            # so, we're typing it as `str` here
            if schema.index is not None:
                index = pandera.Index(
                    schema.index.dtype
                    if not schema.index.dtype.startswith("cat")
                    else str
                )
            else:
                index = None

            self._pandera_schema = pandera.DataFrameSchema(
                pandera_columns,
                coerce=schema.coerce_dtype,
                strict=schema.maximal_set,
                ordered=schema.ordered_set,
                index=index,
            )
        if (
            schema.itype == "Composite"
        ):  # backward compat, should be migrated to Feature.name
            columns_field = Feature.name
        else:
            columns_field = parse_cat_dtype(schema.itype, is_itype=True)["field"]
        # in the DataFrameCatManager, we use the
        # actual columns of the dataset, not the pandera columns
        # the pandera columns might have additional optional columns
        self._cat_manager = DataFrameCatManager(
            self._dataset,
            columns_field=columns_field,
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
        if self._artifact is not None:
            raise RuntimeError(
                "Cannot mutate the dataset when an artifact is passed! Please load the dataset into memory using `dataset.load()` and pass it to a curator."
            )

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
        if self._pandera_schema is not None:
            try:
                # first validate through pandera
                self._pandera_schema.validate(self._dataset, lazy=True)
                # then validate lamindb categoricals
                self._cat_manager_validate()
            except (pandera.errors.SchemaError, pandera.errors.SchemaErrors) as err:
                self._is_validated = False
                has_dtype_error = "WRONG_DATATYPE" in str(err)
                error_msg = str(err)
                if has_dtype_error:
                    error_msg += "   ▶ Hint: Consider setting 'coerce_dtype=True' to attempt coercing/converting values during validation to the pre-defined dtype."
                raise ValidationError(error_msg) from err
        else:
            self._cat_manager_validate()


class DataFrameCurator(SlotsCurator):
    # the example in the docstring is tested in test_curators_quickstart_example
    """Curator for `DataFrame`.

    Args:
        dataset: The DataFrame-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.
        slot: Indicate the slot in a composite curator for a composite data structure.

    Examples:

        For a simple example using a flexible schema, see :meth:`~lamindb.Artifact.from_dataframe`.

        Here is an example that enforces a minimal set of columns in the dataframe.

        .. literalinclude:: scripts/curate_dataframe_minimal_errors.py
            :language: python

        Under-the-hood, this used the following schema.

        .. literalinclude:: scripts/define_mini_immuno_schema_flexible.py
            :language: python

        Valid features & labels were defined as:

        .. literalinclude:: scripts/define_mini_immuno_features_labels.py
            :language: python

        It is also possible to curate the `attrs` slot.

        .. literalinclude:: scripts/curate_dataframe_attrs.py
            :language: python
    """

    def __init__(
        self,
        dataset: pd.DataFrame | Artifact,
        schema: Schema,
        *,
        slot: str | None = None,
        features: dict[str, Any] | None = None,
    ) -> None:
        # loads or opens dataset, dataset may be an artifact
        super().__init__(dataset=dataset, schema=schema, features=features)
        # uses open dataset at self._dataset
        self._atomic_curator = ComponentCurator(
            dataset=self._dataset,
            schema=schema,
            slot=slot,
        )
        # Handle (nested) attrs
        if slot is None and schema.slots:
            for slot_name, slot_schema in schema.slots.items():
                if slot_name.startswith("attrs"):
                    path_parts = slot_name.split(":")
                    attrs_dict = getattr(self._dataset, "attrs", None)
                    if attrs_dict is not None:
                        if len(path_parts) == 1:
                            data = attrs_dict
                        else:
                            deeper_keys = path_parts[1:]
                            data = _resolve_schema_slot_path(
                                attrs_dict, deeper_keys, slot_name, "attrs"
                            )
                        df = convert_dict_to_dataframe_for_validation(data, slot_schema)
                        self._slots[slot_name] = ComponentCurator(
                            df, slot_schema, slot=slot_name
                        )
                elif slot_name != "__external__":
                    raise ValueError(
                        f"Slot '{slot_name}' is not supported for DataFrameCurator. Must be 'attrs'."
                    )

    @property
    def cat(self) -> DataFrameCatManager:
        """Manage categoricals by updating registries."""
        return self._atomic_curator.cat

    def standardize(self) -> None:
        """Standardize the dataset.

        - Adds missing columns for features
        - Fills missing values for features with default values
        """
        self._atomic_curator.standardize()
        for slot_curator in self._slots.values():
            slot_curator.standardize()

    @doc_args(VALIDATE_DOCSTRING)
    def validate(self) -> None:
        """{}."""
        self._atomic_curator.validate()
        self._is_validated = self._atomic_curator._is_validated
        super().validate()

    @doc_args(SAVE_ARTIFACT_DOCSTRING)
    def save_artifact(
        self, *, key=None, description=None, revises=None, run=None
    ) -> Artifact:
        """{}."""
        if not self._is_validated:
            self.validate()
        self._slots["columns"] = self._atomic_curator
        try:
            return super().save_artifact(
                key=key, description=description, revises=revises, run=run
            )
        finally:
            del self._slots["columns"]


class ExperimentalDictCurator(DataFrameCurator):
    """Curator for `dict` based on `DataFrameCurator`."""

    def __init__(
        self,
        dataset: dict | Artifact,
        schema: Schema,
        slot: str | None = None,
    ) -> None:
        if not isinstance(dataset, dict) and not isinstance(dataset, Artifact):
            raise InvalidArgument("The dataset must be a dict or dict-like artifact.")
        if isinstance(dataset, Artifact):
            assert dataset.otype == "dict", "Artifact must be of otype 'dict'."  # noqa: S101
            d = dataset.load(is_run_input=False)
        else:
            d = dataset
        df = convert_dict_to_dataframe_for_validation(d, schema)
        super().__init__(df, schema, slot=slot)


def _resolve_schema_slot_path(
    target_dict: dict[str, Any], slot_keys: Iterable[str], slot: str, base_path: str
) -> Any:
    """Resolve a schema slot path by traversing nested dictionary keys.

    Args:
        target_dict: Root dictionary to traverse
        slot_keys: Sequence of keys defining the paths to traverse
        slot_name: Schema slot identifier for error context
        base_path: Base path string for error context

    Returns:
        The value at the resolved path
    """
    current = target_dict

    for key in slot_keys:
        base_path += f"['{key}']"
        try:
            current = current[key]
        except (
            KeyError,
            TypeError,
        ):  # if not a dict, raises TypeError; if a dict and key not found, raises KeyError
            available = (
                list(current.keys())
                if isinstance(current, dict)
                else "none (not a dict)"
            )
            raise InvalidArgument(
                f"Schema slot '{slot}' requires keys {base_path} but key '{key}' "
                f"not found. Available keys at this level: {available}."
            ) from None

    return current


def _handle_dict_slots(
    dataset: ScverseDataStructures, slot: str
) -> tuple[pd.DataFrame | None, str | None, str | None]:
    """Handle dict-based slot paths (uns/attrs standalone or of modalities) for all ScverseCurators.

    Supports two patterns:
        - Direct dict access: "uns", "attrs", "uns:key1:key2", "attrs:key"
        - Modality dict access: "modality:uns"

    Args:
        dataset: The scverse datastructure object
        slot: The slot path string to parse like 'uns:path:to'.

    Returns:
        tuple: (dataframe, modality_key, remaining_slot_path)
            - dataframe: Single-row DataFrame containing the resolved data
            - modality_key: Modality identifier if slot targets modality dict, else None
            - remaining_slot_path: The dict attribute and nested keys as string
    """
    path_parts = slot.split(":")

    # Handle direct dict slots: "uns", "attrs", "uns:key1:key2:..."
    if len(path_parts) >= 1 and path_parts[0] in ["uns", "attrs"]:
        dict_attr = getattr(dataset, path_parts[0], None)
        if dict_attr is not None:
            if len(path_parts) == 1:
                return pd.DataFrame([dict_attr]), None, path_parts[0]

            deeper_keys = path_parts[1:]
            data = _resolve_schema_slot_path(
                dict_attr, deeper_keys, slot, path_parts[0]
            )
            return pd.DataFrame([data]), None, ":".join(path_parts[1:])

    # Handle modality dict slots: "modality:uns", "modality:uns:key1:key2"
    elif len(path_parts) >= 2 and path_parts[1] in ["uns", "attrs"]:
        modality, dict_name = path_parts[0], path_parts[1]
        try:
            modality_dataset = dataset[modality]
            dict_attr = getattr(modality_dataset, dict_name, None)
            if dict_attr is not None:
                if len(path_parts) == 2:
                    return pd.DataFrame([dict_attr]), modality, dict_name

                deeper_keys = path_parts[2:]
                data = _resolve_schema_slot_path(
                    dict_attr, deeper_keys, slot, f"{modality}.{dict_name}"
                )
                return pd.DataFrame([data]), modality, ":".join(path_parts[1:])
        except (KeyError, AttributeError):
            pass
    else:
        raise InvalidArgument(
            f"Invalid dict slot pattern '{slot}'. Expected formats: "
            f"'uns', 'attrs', 'uns:key', 'attrs:key', 'modality:uns'"
        )

    return None, None, None


@doc_args(SLOTS_DETAILS_DOCSTRING)
class AnnDataCurator(SlotsCurator):
    """Curator for `AnnData`.

    {}

    Args:
        dataset: The AnnData-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Examples:

        Curate Ensembl gene IDs and valid features in obs:

        .. literalinclude:: scripts/curate_anndata_flexible.py
            :language: python
            :caption: curate_anndata_flexible.py

        Curate `uns` dictionary:

        .. literalinclude:: scripts/curate_anndata_uns.py
            :language: python
            :caption: curate_anndata_uns.py

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

        for slot, slot_schema in schema.slots.items():
            if slot not in {"var", "var.T", "obs"} and not slot.startswith("uns"):
                raise ValueError(
                    f"AnnDataCurator currently only supports the slots 'var', 'var.T', 'obs', and 'uns', not {slot}"
                )
            if slot.startswith("uns"):
                df, _, _ = _handle_dict_slots(self._dataset, slot)
            elif slot in {"obs", "var", "var.T"}:
                df = (
                    getattr(self._dataset, slot.strip(".T")).T
                    if slot == "var.T"
                    or (
                        slot == "var"
                        and schema.slots["var"].itype not in {None, "Feature"}
                    )
                    else getattr(self._dataset, slot)
                )
            self._slots[slot] = ComponentCurator(df, slot_schema, slot=slot)

            # Handle var index naming for backward compat
            if slot == "var" and schema.slots["var"].itype not in {None, "Feature"}:
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
    slots: dict[str, ComponentCurator],
) -> None:
    """Assigns var_fields and categoricals for multimodal data curators."""
    if modality is not None:
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


@doc_args(SLOTS_DETAILS_DOCSTRING)
class MuDataCurator(SlotsCurator):
    """Curator for `MuData`.

    {}

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
            # Handle slots: "mdata.uns", "modality:uns"
            if "uns" in slot:
                df, modality, modality_slot = _handle_dict_slots(self._dataset, slot)
            else:
                # Handle slots: "modality:obs", "modality:var"
                parts = slot.split(":")
                if len(parts) == 2:
                    modality, modality_slot = parts
                    try:
                        schema_dataset = self._dataset[modality]
                        df = getattr(schema_dataset, modality_slot.rstrip(".T"))
                    except KeyError:
                        raise InvalidArgument(
                            f"Modality '{modality}' not found in MuData"
                        ) from None
                    except AttributeError:
                        raise InvalidArgument(
                            f"Attribute '{modality_slot}' not found on modality '{modality}'"
                        ) from None
                else:
                    # Handle slots: "mdata:obs", "mdata:var" (uns is a dictionary and gets handled above)
                    modality, modality_slot = None, slot
                    schema_dataset = self._dataset
                    df = getattr(schema_dataset, modality_slot.rstrip(".T"))

            # Transpose var if necessary
            if modality_slot == "var" and schema.slots[slot].itype not in {
                None,
                "Feature",
            }:
                logger.warning(
                    "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                )
                df = df.T
            elif modality_slot == "var.T":
                df = df.T

            self._slots[slot] = ComponentCurator(df, slot_schema, slot=slot)

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


@doc_args(SLOTS_DETAILS_DOCSTRING)
class SpatialDataCurator(SlotsCurator):
    """Curator for `SpatialData`.

    {}

    Args:
        dataset: The SpatialData-like object to validate & annotate.
        schema: A :class:`~lamindb.Schema` object that defines the validation constraints.

    Example:
        .. literalinclude:: scripts/curate_spatialdata.py
            :language: python
            :caption: curate_spatialdata.py

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
            # Handle slots: "sdata:attrs"
            if slot.startswith("attrs"):
                df, table_key, table_slot = _handle_dict_slots(self._dataset, slot)
            else:
                parts = slot.split(":")
                # Handle slots: "tables:table_key:obs", "tables:table_key:var"
                if len(parts) == 3 and parts[0] == "tables":
                    table_key, table_slot = parts[1], parts[2]
                    try:
                        slot_object = self._dataset.tables[table_key]
                        df = getattr(slot_object, table_slot.rstrip(".T"))
                    except KeyError:
                        raise InvalidArgument(
                            f"Table '{table_key}' not found in sdata.tables"
                        ) from None
                    except AttributeError:
                        raise InvalidArgument(
                            f"Attribute '{table_slot}' not found on table '{table_key}'"
                        ) from None
                else:
                    # Handle legacy single keys for backward compatibility
                    if len(parts) == 1 and parts[0] != "attrs":
                        logger.warning(
                            f"please prefix slot {slot} with 'attrs:' going forward"
                        )
                        try:
                            df = pd.DataFrame([self._dataset.attrs[slot]])
                            table_key = None
                            table_slot = slot
                        except KeyError:
                            raise InvalidArgument(
                                f"Slot '{slot}' not found in sdata.attrs"
                            ) from None
                    else:
                        raise InvalidArgument(f"Unrecognized slot format: {slot}")

            # Handle var transposition logic
            if table_slot == "var" and schema.slots[slot].itype not in {
                None,
                "Feature",
            }:
                logger.warning(
                    "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                )
                df = df.T
            elif table_slot == "var.T":
                df = df.T

            self._slots[slot] = ComponentCurator(df, slot_schema, slot)

            _assign_var_fields_categoricals_multimodal(
                modality=table_key,
                slot_type=table_slot,
                slot=slot,
                slot_schema=slot_schema,
                var_fields=self._var_fields,
                cat_vectors=self._cat_vectors,
                slots=self._slots,
            )

        self._columns_field = self._var_fields


@doc_args(SLOTS_DETAILS_DOCSTRING)
class TiledbsomaExperimentCurator(SlotsCurator):
    """Curator for `tiledbsoma.Experiment`.

    {}

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
                _, modality_slot = slot.split(":")
                schema_dataset = (
                    self._dataset.ms[modality_slot.removesuffix(".T")]
                    .var.read()
                    .concat()
                    .to_pandas()
                    .drop("soma_joinid", axis=1, errors="ignore")
                )

                self._slots[slot] = ComponentCurator(
                    (schema_dataset.T if modality_slot == "var.T" else schema_dataset),
                    slot_schema,
                )
            else:
                # global Experiment obs slot
                modality_slot = slot
                schema_dataset = (
                    self._dataset.obs.read()
                    .concat()
                    .to_pandas()
                    .drop(["soma_joinid", "obs_id"], axis=1, errors="ignore")
                )
                self._slots[slot] = ComponentCurator(
                    schema_dataset,
                    slot_schema,
                )

            _assign_var_fields_categoricals_multimodal(
                modality=slot,  # not passing `measurement` here because it's a constant. The slot has the actual modality
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
        filter_str: str = "",
        subtypes_list: list[str] = None,
        maximal_set: bool = True,  # whether unvalidated categoricals cause validation failure.
    ) -> None:
        if subtypes_list is None:
            subtypes_list = []
        self._values_getter = values_getter
        self._values_setter = values_setter
        self._field = field
        self._key = key
        self._source = source
        self._organism = None
        self._validated: None | list[str] = None
        self._non_validated: None | list[str] = None
        self._synonyms: None | dict[str, str] = None
        self._filter_str = filter_str
        self._subtypes_list = subtypes_list
        self._subtype_query_set = None
        self._cat_manager = cat_manager
        self.feature = feature
        self.records = None
        self._maximal_set = maximal_set
        self._type_record = None

        self._all_filters = {"source": self._source, "organism": self._organism}

        if self._filter_str:
            self._all_filters.update(
                resolve_relation_filters(
                    parse_filter_string(self._filter_str), self._field.field.model
                )  # type: ignore
            )

        # get the dtype associated record based on the nested subtypes
        if self._subtypes_list:
            type_filters = {"name": self._subtypes_list[-1]}
            if len(self._subtypes_list) > 1:
                for i, nested_subtype in enumerate(reversed(self._subtypes_list[:-1])):
                    filter_key = f"{'type__' * (i + 1)}name"
                    type_filters[filter_key] = nested_subtype
            try:
                self._type_record = self._field.field.model.get(**type_filters)
            except Exception as e:
                raise InvalidArgument(
                    f"Error retrieving type record with filters {type_filters} for field {self._field.field.name}."
                ) from e

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
            logger.warning(f"no values were validated for {self._key}!")
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
        filter_kwargs = get_current_filter_kwargs(registry, self._all_filters)

        valid_from_values_kwargs = {}
        for key, value in filter_kwargs.items():
            if key in {"field", "organism", "source", "mute"}:
                valid_from_values_kwargs[key] = value
            elif hasattr(registry, key) and "__" not in key:
                valid_from_values_kwargs[key] = value

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

        # if values are SQLRecord, we don't need to validate them
        if all(isinstance(v, SQLRecord) for v in str_values):
            assert all(v._state.adding is False for v in str_values), (
                "All records must be saved."
            )
            self.records = str_values  # type: ignore
            validated_labels = str_values  # type: ignore
            return validated_labels, []

        # inspect the default instance and save validated records from public
        if self._type_record is not None:
            related_name = registry._meta.get_field("type").remote_field.related_name
            if registry.__name__ == "Record":
                self._subtype_query_set = self._type_record.query_records()
            else:
                self._subtype_query_set = getattr(self._type_record, related_name).all()
            values_array = np.array(str_values)
            validated_mask = self._subtype_query_set.validate(  # type: ignore
                values_array, field=self._field, **filter_kwargs, mute=True
            )
            validated_labels, non_validated_labels = (
                values_array[validated_mask],
                values_array[~validated_mask],
            )
            records = registry.from_values(
                validated_labels,
                field=self._field,
                **valid_from_values_kwargs,
                mute=True,
            )
        else:
            existing_and_public_records = registry.from_values(
                str_values, field=self._field, **valid_from_values_kwargs, mute=True
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
            non_validated_records = Feature.from_dataframe(df.loc[:, nonval_columns])
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
                if self._type_record is not None:
                    # if type_record is set, we need to set the type for new records
                    init_kwargs["type"] = self._type_record
                # here we create non-validated records skipping validation since we already ensured that they don't exist
                non_validated_records.append(
                    registry(**init_kwargs, **create_kwargs, _skip_validation=True)
                )
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

        kwargs_current = get_current_filter_kwargs(registry, self._all_filters)

        valid_inspect_kwargs = {}
        for key, value in kwargs_current.items():
            if key in {"field", "organism", "source", "mute", "from_source"}:
                valid_inspect_kwargs[key] = value
            elif hasattr(registry, key) and "__" not in key:
                valid_inspect_kwargs[key] = value

        # inspect values from the default instance, excluding public
        registry_or_queryset = registry
        if self._subtype_query_set is not None:
            registry_or_queryset = self._subtype_query_set
        inspect_result = registry_or_queryset.inspect(
            values,
            field=self._field,
            mute=True,
            from_source=False,
            **valid_inspect_kwargs,
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
                **valid_inspect_kwargs,
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
                check_organism = ""
                if registry.__base__.__name__ == "BioRecord":
                    import bionty as bt
                    from bionty._organism import is_organism_required

                    if is_organism_required(registry):
                        organism = (
                            valid_inspect_kwargs.get("organism", False)
                            or bt.settings.organism.name
                        )
                        check_organism = f"fix organism '{organism}', "
                warning_message += f"    → {check_organism}fix typos, remove non-existent values, or save terms via: {colors.cyan(non_validated_hint_print)}"
                if self._subtype_query_set is not None:
                    warning_message += f"\n    → a valid label for subtype '{self._subtypes_list[-1]}' has to be one of {self._subtype_query_set.to_list('name')}"
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

        self._cat_vectors["columns"] = CatVector(
            values_getter=lambda: self._dataset.keys(),  # lambda ensures the inplace update
            values_setter=lambda new_values: setattr(
                self._dataset, "columns", pd.Index(new_values)
            )
            if isinstance(self._dataset, pd.DataFrame)
            else None,
            field=columns_field,
            key="columns" if isinstance(self._dataset, pd.DataFrame) else "keys",
            source=self._sources.get("columns"),
            cat_manager=self,
            maximal_set=self._maximal_set,
        )
        for feature in self._categoricals:
            result = parse_dtype(feature.dtype)[
                0
            ]  # TODO: support composite dtypes for categoricals
            key = feature.name
            self._cat_vectors[key] = CatVector(
                values_getter=lambda k=key: self._dataset[
                    k
                ],  # Capture key as default argument
                values_setter=lambda new_values, k=key: self._dataset.__setitem__(
                    k, new_values
                ),
                field=result["field"],
                key=key,
                source=self._sources.get(key),
                feature=feature,
                cat_manager=self,
                filter_str=result["filter_str"],
                subtypes_list=result["subtypes_list"],
            )
        if index is not None and index.dtype.startswith("cat"):
            result = parse_dtype(index.dtype)[0]
            key = "index"
            self._cat_vectors[key] = CatVector(
                values_getter=self._dataset.index,
                values_setter=lambda new_values: setattr(
                    self._dataset, "index", new_values
                ),
                field=result["field"],
                key=key,
                feature=index,
                cat_manager=self,
                filter_str=result["filter_str"],
                subtypes_list=result["subtypes_list"],
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

    def __repr__(self) -> str:
        cls_name = colors.green(self.__class__.__name__)

        status_str = (
            f"{colors.green('validated')}"
            if self._is_validated
            else f"{colors.yellow('unvalidated')}"
        )

        info_parts = []

        cat_count = len(self._categoricals)
        if cat_count > 0:
            info_parts.append(f"categorical_features={cat_count}")

        if self._slot:
            info_parts.append(f"slot: {colors.italic(self._slot)}")

        info_str = ", ".join(info_parts)
        if info_str:
            return f"{cls_name}({info_str}, {status_str})"
        else:
            return f"{cls_name}({status_str})"

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
            raise RuntimeError(
                "Cannot mutate the dataset when an artifact is passed! Please load the dataset into memory using `dataset.load()` and pass it to a curator."
            )

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


def get_current_filter_kwargs(
    registry: type[SQLRecord], kwargs: dict[str, SQLRecord]
) -> dict:
    """Make sure the source and organism are saved in the same database as the registry."""
    db = registry.filter().db
    filter_kwargs = kwargs.copy()

    for key, value in kwargs.items():
        if isinstance(value, SQLRecord) and value._state.db != "default":
            if db is None or db == "default":
                value_default = copy.copy(value)
                value_default.save()
                filter_kwargs[key] = value_default

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
    curator: SlotsCurator | None = None,
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
            or cat_vector.records is None
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
    if (
        artifact.otype == "DataFrame" and getattr(curator, "_schema", None) is None
    ):  # Prevent overwriting user-defined schemas that contain slots
        features = cat_vectors["columns"].records
        if features is not None:
            index_feature = artifact.schema.index
            feature_set = Schema(
                features=[f for f in features if f != index_feature],
                itype=artifact.schema.itype,
                index=index_feature,
                minimal_set=artifact.schema.minimal_set,
                maximal_set=artifact.schema.maximal_set,
                coerce_dtype=artifact.schema.coerce_dtype,
                ordered_set=artifact.schema.ordered_set,
            )
            if (
                feature_set._state.adding
                and len(features) > settings.annotation.n_max_records
            ):
                logger.important(
                    f"not annotating with {len(features)} features as it exceeds {settings.annotation.n_max_records} (ln.settings.annotation.n_max_records)"
                )
                itype = (
                    Feature.name
                    if artifact.schema.itype == "Composite"  # backward compat
                    else parse_cat_dtype(artifact.schema.itype, is_itype=True)["field"]
                )
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
            validating_schema = slot_curator._schema
            index_feature = validating_schema.index
            feature_set = Schema(
                features=[f for f in features if f != index_feature],
                itype=validating_schema.itype,
                index=index_feature,
                minimal_set=validating_schema.minimal_set,
                maximal_set=validating_schema.maximal_set,
                coerce_dtype=validating_schema.coerce_dtype,
                ordered_set=validating_schema.ordered_set,
            )
            if (
                feature_set._state.adding
                and len(features) > settings.annotation.n_max_records
            ):
                logger.important(
                    f"not annotating with {len(features)} features for slot {slot} as it exceeds {settings.annotation.n_max_records} (ln.settings.annotation.n_max_records)"
                )
                itype = (
                    Feature.name
                    if artifact.schema.slots[slot].itype
                    == "Composite"  # backward compat
                    else parse_cat_dtype(
                        artifact.schema.slots[slot].itype, is_itype=True
                    )["field"]
                )
                feature_set = Schema(itype=itype, n=len(features))
            artifact.feature_sets.add(
                feature_set.save(), through_defaults={"slot": slot}
            )

    slug = ln_setup.settings.instance.slug
    if ln_setup.settings.instance.is_remote:  # pdagma: no cover
        ui_url = ln_setup.settings.instance.ui_url
        logger.important(f"go to {ui_url}/{slug}/artifact/{artifact.uid}")
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
