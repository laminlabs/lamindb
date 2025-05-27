from __future__ import annotations

from typing import TYPE_CHECKING, Any, Type, overload

import numpy as np
from django.db import models
from django.db.models import CASCADE, PROTECT, ManyToManyField
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_string
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

from lamindb.base import ids
from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    IntegerField,
    JSONField,
)
from lamindb.base.types import FieldAttr, ListLike
from lamindb.errors import FieldValidationError, InvalidArgument
from lamindb.models.feature import parse_cat_dtype

from ..errors import ValidationError
from ._describe import format_rich_tree, highlight_time
from ._relations import (
    dict_related_model_to_related_name,
    get_related_name,
)
from .can_curate import CanCurate
from .feature import (
    Feature,
    serialize_dtype,
    serialize_pandas_dtype,
)
from .run import TracksRun, TracksUpdates
from .sqlrecord import (
    BaseSQLRecord,
    IsLink,
    Registry,
    SQLRecord,
    _get_record_kwargs,
    init_self_from_db,
    update_attributes,
)

if TYPE_CHECKING:
    import pandas as pd
    from django.db.models.query_utils import DeferredAttribute

    from .artifact import Artifact
    from .project import Project
    from .query_set import QuerySet


NUMBER_TYPE = "num"
DICT_KEYS_TYPE = type({}.keys())  # type: ignore


def validate_features(features: list[SQLRecord]) -> SQLRecord:
    """Validate and return feature type."""
    try:
        if len(features) == 0:
            raise ValueError("Provide list of features with at least one element")
    except TypeError:
        raise ValueError(
            "Please pass a ListLike of features, not a single feature"
        ) from None
    if not hasattr(features, "__getitem__"):
        raise TypeError("features has to be list-like")
    if not isinstance(features[0], SQLRecord):
        raise TypeError(
            "features has to store feature records! use .from_values() otherwise"
        )
    feature_types = {feature.__class__ for feature in features}
    if len(feature_types) > 1:
        raise TypeError("schema can only contain a single type")
    for feature in features:
        if feature._state.adding:
            raise ValueError("Can only construct feature sets from validated features")
    return next(iter(feature_types))  # return value in set of cardinality 1


def get_features_config(
    features: list[SQLRecord] | tuple[SQLRecord, dict],
) -> tuple[list[SQLRecord], list[tuple[SQLRecord, dict]]]:
    """Get features and their config from the return of feature.with_config()."""
    features_list = []
    configs = []
    try:
        for feature in features:
            if isinstance(feature, tuple):
                features_list.append(feature[0])
                configs.append(feature)  # store the tuple in configs
            else:
                features_list.append(feature)
        return features_list, configs  # type: ignore
    except TypeError:
        return features, configs  # type: ignore


def describe_schema(self: Schema) -> Tree:
    """Create a rich tree visualization of a Schema with its features."""
    otype = self.otype if hasattr(self, "otype") and self.otype else ""
    tree = Tree(
        Text.assemble((self.__class__.__name__, "bold"), (f" {otype}", "bold dim")),
        guide_style="dim",  # dim the connecting lines
    )

    tree.add(f".uid = '{self.uid}'")
    tree.add(f".name = '{self.name}'")
    if self.description:
        tree.add(f".description = '{self.description}'")
    if self.itype:
        tree.add(f".itype = '{self.itype}'")
    if self.type:
        tree.add(f".type = '{self.type}'")
    tree.add(f".ordered_set = {self.ordered_set}")
    tree.add(f".maximal_set = {self.maximal_set}")
    if hasattr(self, "created_by") and self.created_by:
        tree.add(
            Text.assemble(
                ".created_by = ",
                (
                    self.created_by.handle
                    if self.created_by.name is None
                    else f"{self.created_by.handle} ({self.created_by.name})"
                ),
            )
        )
    if hasattr(self, "created_at") and self.created_at:
        tree.add(Text.assemble(".created_at = ", highlight_time(str(self.created_at))))

    members = self.members

    # Add features section
    features = tree.add(
        Text.assemble(
            (self.itype, "violet"),
            (" • ", "dim"),
            (str(members.count()), "pink1"),
        )
    )

    if hasattr(self, "members") and self.members.count() > 0:
        # create a table for the features
        feature_table = Table(
            show_header=True, header_style="dim", box=None, pad_edge=False
        )

        # Add columns
        feature_table.add_column("name", style="", no_wrap=True)
        feature_table.add_column("dtype", style="", no_wrap=True)
        feature_table.add_column("optional", style="", no_wrap=True)
        feature_table.add_column("nullable", style="", no_wrap=True)
        feature_table.add_column("coerce_dtype", style="", no_wrap=True)
        feature_table.add_column("default_value", style="", no_wrap=True)

        # Add rows for each member
        optionals = self.optionals.get()
        for member in self.members:
            feature_table.add_row(
                member.name,
                Text(
                    str(member.dtype)
                ),  # needs to be wrapped in Text to display correctly
                "✓" if optionals.filter(uid=member.uid).exists() else "✗",
                "✓" if member.nullable else "✗",
                "✓" if member.coerce_dtype else "✗",
                str(member.default_value) if member.default_value else "unset",
            )

        # Add the table to the features branch
        features.add(feature_table)

    return tree


class SchemaOptionals:
    """Manage and access optional features in a schema."""

    def __init__(self, schema) -> None:
        self.schema = schema

    def get_uids(self) -> list[str]:
        """Get the uids of the optional features.

        Does **not** need an additional query to the database, while `get()` does.
        """
        if (
            self.schema._aux is not None
            and "af" in self.schema._aux
            and "1" in self.schema._aux["af"]
        ):
            return self.schema._aux["af"]["1"]
        else:
            return []

    def get(self) -> QuerySet:
        """Get the optional features."""
        uids = self.get_uids()
        if uids:
            return Feature.objects.filter(uid__in=uids).order_by("links_schema__id")
        else:
            return Feature.objects.none()  # empty QuerySet

    def set(self, features: list[Feature]) -> None:
        """Set the optional features (overwrites whichever schemas are currently optional)."""
        if not isinstance(features, list) or not all(
            isinstance(f, Feature) for f in features
        ):
            raise TypeError("features must be a list of Feature records!")
        self.schema._aux = self.schema._aux or {}
        if len(features) > 0:
            self.schema._aux.setdefault("af", {})["1"] = [f.uid for f in features]

    def remove(self, features: Feature | list[Feature]) -> None:
        """Make one or multiple features required by removing them from the set of optional features."""
        if not isinstance(features, list):
            features = [features]
        if not all(isinstance(f, Feature) for f in features):
            raise TypeError("features must be a list of Feature records!")
        if len(features) > 0:
            self.schema._aux = self.schema._aux or {}
            if "1" in self.schema._aux.get("af", {}):
                for feature in features:
                    self.schema._aux["af"]["1"].remove(feature.uid)

    def add(self, features: Feature | list[Feature]) -> None:
        """Make one or multiple features optional by adding them to the set of optional features."""
        self.schema._aux = self.schema._aux or {}
        if not isinstance(features, list):
            features = [features]
        if not all(isinstance(f, Feature) for f in features):
            raise TypeError("features must be a list of Feature records!")
        if len(features) > 0:
            if "1" not in self.schema._aux.setdefault("af", {}):
                self.set(features)
            else:
                self.schema._aux.setdefault("af", {})["1"].extend(
                    [f.uid for f in features]
                )


KNOWN_SCHEMAS = {
    "kMi7B_N88uu-YnbTLDU-DA": "0000000000000000",  # valid_features
    "1gocc_TJ1RU2bMwDRK-WUA": "0000000000000001",  # valid_ensembl_gene_ids
    "GTxxM36n9tocphLfdbNt9g": "0000000000000002",  # anndata_ensembl_gene_ids_and_valid_features_in_obs
}


class Schema(SQLRecord, CanCurate, TracksRun):
    """Schemas of a dataset such as the set of columns of a `DataFrame`.

    Composite schemas can have multiple slots, e.g., for an `AnnData`, one schema for slot `obs` and another one for `var`.

    Args:
        features: `list[SQLRecord] | list[tuple[Feature, dict]] | None = None` Feature
            records, e.g., `[Feature(...), Feature(...)]` or Features with their config, e.g., `[Feature(...).with_config(optional=True)]`.
        index: `Feature | None = None` A :class:`~lamindb.Feature` record to validate an index of a `DataFrame` and therefore also, e.g., `AnnData` obs and var indices.
        slots: `dict[str, Schema] | None = None` A dictionary mapping slot names to :class:`~lamindb.Schema` objects.
        name: `str | None = None` Name of the Schema.
        description: `str | None = None` Description of the Schema.
        flexible: `bool | None = None` Whether to include any feature of the same `itype` in validation
            and annotation. If no Features are passed, defaults to `True`, otherwise to `False`.
            This means that if you explicitly pass Features, any additional Features will be disregarded during validation & annotation.
        type: `Schema | None = None` Type of Schema to group measurements by.
            Define types like `ln.Schema(name="ProteinPanel", is_type=True)`.
        is_type: `bool = False` Whether the Schema is a Type.
        itype: `str | None = None` The feature identifier type (e.g. :class:`~lamindb.Feature`, :class:`~bionty.Gene`, ...).
        otype: `str | None = None` An object type to define the structure of a composite schema (e.g., DataFrame, AnnData).
        dtype: `str | None = None` The simple type (e.g., "num", "float", "int").
            Defaults to `None` for sets of :class:`~lamindb.Feature` records and to `"num"` (e.g., for sets of :class:`~bionty.Gene`) otherwise.
        minimal_set: `bool = True` Whether all passed Features are required by default.
            See :attr:`~lamindb.Schema.optionals` for more-fine-grained control.
        maximal_set: `bool = False` Whether additional Features are allowed.
        ordered_set: `bool = False` Whether Features are required to be ordered.
        coerce_dtype: `bool = False` When True, attempts to coerce values to the specified dtype
            during validation, see :attr:`~lamindb.Schema.coerce_dtype`.

    See Also:
        :meth:`~lamindb.Artifact.from_df`
            Validate & annotate a `DataFrame` with a schema.
        :meth:`~lamindb.Artifact.from_anndata`
            Validate & annotate an `AnnData` with a schema.
        :meth:`~lamindb.Artifact.from_mudata`
            Validate & annotate an `MuData` with a schema.
        :meth:`~lamindb.Artifact.from_spatialdata`
            Validate & annotate a `SpatialData` with a schema.

    Examples:

        The typical way to create a schema::

            import lamindb as ln
            import bionty as bt
            import pandas as pd

            # a schema with a single required feature
            schema = ln.Schema(
                features=[
                    ln.Feature(name="required_feature", dtype=str).save(),
                ],
            ).save()

            # a schema that constrains feature identifiers to be a valid ensembl gene ids or feature names
            schema = ln.Schema(itype=bt.Gene.ensembl_gene_id)
            schema = ln.Schema(itype=ln.Feature)  # is equivalent to itype=ln.Feature.name

            # a schema that requires a single feature but also validates & annotates any additional features with valid feature names
            schema = ln.Schema(
                features=[
                    ln.Feature(name="required_feature", dtype=str).save(),
                ],
                itype=ln.Schema(itype=ln.Feature),
                flexible=True,
            ).save()

        Passing options to the `Schema` constructor::

            # also validate the index
            schema = ln.Schema(
                features=[
                    ln.Feature(name="required_feature", dtype=str).save(),
                ],
                index=ln.Feature(name="sample", dtype=ln.ULabel).save(),
            ).save()

            # mark a single feature as optional and ignore other features of the same identifier type
            schema = ln.Schema(
                features=[
                    ln.Feature(name="required_feature", dtype=str).save(),
                    ln.Feature(name="feature2", dtype=int).save().with_config(optional=True),
                ],
            ).save()

        Alternative constructors (:meth:`~lamindb.Schema.from_values`, :meth:`~lamindb.Schema.from_df`)::

            # parse & validate identifier values
            schema = ln.Schema.from_values(
                adata.var["ensemble_id"],
                field=bt.Gene.ensembl_gene_id,
                organism="mouse",
            ).save()

            # from a dataframe
            df = pd.DataFrame({"feat1": [1, 2], "feat2": [3.1, 4.2], "feat3": ["cond1", "cond2"]})
            schema = ln.Schema.from_df(df)
    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"
    _aux_fields: dict[str, tuple[str, type]] = {
        "0": ("coerce_dtype", bool),
        "1": ("optionals", list[str]),
        "2": ("flexible", bool),
        "3": ("index_feature_uid", str),
    }

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    # Before lamindb 1.5, it was 20 char long. Since lamindb 1.5, it is 16 char long.
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=20)
    """A universal id."""
    name: str | None = CharField(max_length=150, null=True, db_index=True)
    """A name."""
    description: str | None = CharField(null=True, db_index=True)
    """A description."""
    n: int = IntegerField()
    """Number of features in the schema."""
    type: Schema | None = ForeignKey("self", PROTECT, null=True, related_name="schemas")
    """Type of schema.

    Allows to group schemas by type, e.g., all meassurements evaluating gene expression vs. protein expression vs. multi modal.

    You can define types via `ln.Schema(name="ProteinPanel", is_type=True)`.

    Here are a few more examples for type names: `'ExpressionPanel'`, `'ProteinPanel'`, `'Multimodal'`, `'Metadata'`, `'Embedding'`.
    """
    instances: Schema
    """Schemas of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    itype: str | None = CharField(
        max_length=120, db_index=True, null=True, editable=False
    )
    """A registry that stores feature identifier types used in this schema, e.g., `'Feature'` or `'bionty.Gene'`.

    Depending on `itype`, `.members` stores, e.g., `Feature` or `bionty.Gene` records.
    """
    otype: str | None = CharField(max_length=64, db_index=True, null=True)
    """Default Python object type, e.g., DataFrame, AnnData."""
    dtype: str | None = CharField(max_length=64, null=True, editable=False)
    """Data type, e.g., "num", "float", "int". Is `None` for :class:`~lamindb.Feature`.

    For :class:`~lamindb.Feature`, types are expected to be heterogeneous and defined on a per-feature level.
    """
    hash: str | None = CharField(
        max_length=HASH_LENGTH, db_index=True, null=True, editable=False
    )
    """A hash of the set of feature identifiers.

    For a composite schema, the hash of hashes.
    """
    minimal_set: bool = BooleanField(default=True, db_index=True, editable=False)
    """Whether all passed features are to be considered required by default (default `True`).

    Note that features that are explicitly marked as `optional` via `feature.with_config(optional=True)`
    are **not** required even if this `minimal_set` is true.
    """
    ordered_set: bool = BooleanField(default=False, db_index=True, editable=False)
    """Whether features are required to be ordered (default `False`)."""
    maximal_set: bool = BooleanField(default=False, db_index=True, editable=False)
    """Whether all features present in the dataset must be in the schema (default `False`).

    If `False`, additional features are allowed to be present in the dataset.

    If `True`, no additional features are allowed to be present in the dataset.
    """
    components: Schema = ManyToManyField(
        "self", through="SchemaComponent", symmetrical=False, related_name="composites"
    )
    """Components of this schema."""
    composites: Schema
    """The composite schemas that contains this schema as a component.

    For example, an `AnnData` composes multiple schemas: `var[DataFrameT]`, `obs[DataFrame]`, `obsm[Array]`, `uns[dict]`, etc.
    """
    features: Feature
    """The features contained in the schema."""
    artifacts: Artifact
    """The artifacts that measure a feature set that matches this schema."""
    validated_artifacts: Artifact
    """The artifacts that were validated against this schema with a :class:`~lamindb.curators.core.Curator`."""
    projects: Project
    """Linked projects."""
    _curation: dict[str, Any] = JSONField(default=None, db_default=None, null=True)
    # lamindb v2
    # _itype: ContentType = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    # ""Index of the registry that stores the feature identifiers, e.g., `Feature` or `Gene`."""
    # -- the following two fields are dynamically removed from the API for now
    validated_by: Schema | None = ForeignKey(
        "self", PROTECT, related_name="validated_schemas", default=None, null=True
    )
    # """The schema that validated this schema during curation.

    # When performing validation, the schema that enforced validation is often less concrete than what is validated.

    # For instance, the set of measured features might be a superset of the minimally required set of features.
    # """
    # validated_schemas: Schema
    # """The schemas that were validated against this schema with a :class:`~lamindb.curators.core.Curator`."""
    composite: Schema | None = ForeignKey(
        "self", PROTECT, related_name="+", default=None, null=True
    )
    # The legacy foreign key
    slot: str | None = CharField(max_length=100, db_index=True, null=True)
    # The legacy slot

    @overload
    def __init__(
        self,
        features: list[SQLRecord] | list[tuple[Feature, dict]] | None = None,
        index: Feature | None = None,
        slots: dict[str, Schema] | None = None,
        name: str | None = None,
        description: str | None = None,
        itype: str | Registry | FieldAttr | None = None,
        flexible: bool | None = None,
        type: Schema | None = None,
        is_type: bool = False,
        otype: str | None = None,
        dtype: str | Type[int | float | str] | None = None,  # noqa
        ordered_set: bool = False,
        minimal_set: bool = True,
        maximal_set: bool = False,
        coerce_dtype: bool = False,
        n: int | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None
        if len(args) > 1:
            raise ValueError("Only one non-keyword arg allowed: features")

        features: list[SQLRecord] | None = (
            args[0] if args else kwargs.pop("features", [])
        )
        index: Feature | None = kwargs.pop("index", None)
        slots: dict[str, Schema] = kwargs.pop("slots", {})
        name: str | None = kwargs.pop("name", None)
        description: str | None = kwargs.pop("description", None)
        itype: str | SQLRecord | DeferredAttribute | None = kwargs.pop("itype", None)
        flexible: bool | None = kwargs.pop("flexible", None)
        type: Feature | None = kwargs.pop("type", None)
        is_type: bool = kwargs.pop("is_type", False)
        otype: str | None = kwargs.pop("otype", None)
        dtype: str | None = kwargs.pop("dtype", None)
        minimal_set: bool = kwargs.pop("minimal_set", True)
        ordered_set: bool = kwargs.pop("ordered_set", False)
        maximal_set: bool = kwargs.pop("maximal_set", False)
        coerce_dtype: bool | None = kwargs.pop("coerce_dtype", False)
        using: bool | None = kwargs.pop("using", None)
        n_features: int | None = kwargs.pop("n", None)
        # backward compat
        if not slots:
            if "components" in kwargs:
                logger.warning(
                    "`components` as a keyword argument is deprecated, please use `slots` instead"
                )
                slots = kwargs.pop("components")
        if kwargs:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Schema)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        (
            features,
            validated_kwargs,
            optional_features,
            features_registry,
            flexible,
        ) = self._validate_kwargs_calculate_hash(
            features=features,
            index=index,
            slots=slots,
            name=name,
            description=description,
            itype=itype,
            flexible=flexible,
            type=type,
            is_type=is_type,
            otype=otype,
            dtype=dtype,
            minimal_set=minimal_set,
            ordered_set=ordered_set,
            maximal_set=maximal_set,
            coerce_dtype=coerce_dtype,
            n_features=n_features,
        )
        schema = (
            Schema.objects.using(using)
            .filter(hash=validated_kwargs["hash"])
            .one_or_none()
        )
        if schema is not None:
            logger.important(f"returning existing schema with same hash: {schema}")
            init_self_from_db(self, schema)
            update_attributes(self, validated_kwargs)
            self.optionals.set(optional_features)
            return None
        self._slots: dict[str, Schema] = {}
        if features:
            self._features = (get_related_name(features_registry), features)  # type: ignore
        elif slots:
            for slot_key, component in slots.items():
                if component._state.adding:
                    raise InvalidArgument(
                        f"schema for {slot_key} {component} must be saved before use"
                    )
            self._slots = slots
        if validated_kwargs["hash"] in KNOWN_SCHEMAS:
            validated_kwargs["uid"] = KNOWN_SCHEMAS[validated_kwargs["hash"]]
        else:
            validated_kwargs["uid"] = ids.base62_16()
        super().__init__(**validated_kwargs)
        # manipulating aux fields is easier after calling super().__init__()
        self.optionals.set(optional_features)
        self.flexible = flexible
        if index is not None:
            self._index_feature_uid = index.uid

    def _validate_kwargs_calculate_hash(
        self,
        features: list[SQLRecord],
        index: Feature | None,
        slots: dict[str, Schema],
        name: str | None,
        description: str | None,
        itype: str | SQLRecord | DeferredAttribute | None,
        flexible: bool | None,
        type: Feature | None,
        is_type: bool,
        otype: str | None,
        dtype: str | None,
        minimal_set: bool,
        ordered_set: bool,
        maximal_set: bool,
        coerce_dtype: bool,
        n_features: int | None,
        optional_features_manual: list[Feature] | None = None,
    ) -> tuple[list[Feature], dict[str, Any], list[Feature], Registry, bool]:
        optional_features = []
        features_registry: Registry = None
        if itype is not None:
            if itype != "Composite":
                itype = serialize_dtype(itype, is_itype=True)
        if index is not None:
            if not isinstance(index, Feature):
                raise TypeError("index must be a Feature")
            features.insert(0, index)
        if features:
            features, configs = get_features_config(features)
            features_registry = validate_features(features)
            itype_compare = features_registry.__get_name_with_module__()
            if itype is not None:
                assert itype.startswith(itype_compare), str(itype_compare)  # noqa: S101
            else:
                itype = itype_compare
            if n_features is not None:
                if n_features != len(features):
                    logger.important(f"updating to n {len(features)} features")
            n_features = len(features)
            if features_registry == Feature:
                optional_features = [
                    config[0] for config in configs if config[1].get("optional")
                ]
                if optional_features:
                    assert optional_features_manual is None  # noqa: S101
                if not optional_features and optional_features_manual is not None:
                    optional_features = optional_features_manual
        elif n_features is None:
            n_features = -1
        if dtype is None:
            dtype = None if itype is not None and itype == "Feature" else NUMBER_TYPE
        else:
            dtype = get_type_str(dtype)
        flexible_default = n_features < 0
        if flexible is None:
            flexible = flexible_default
        if slots:
            itype = "Composite"
            if otype is None:
                raise InvalidArgument("Please pass otype != None for composite schemas")
        if itype is not None and not isinstance(itype, str):
            itype_str = serialize_dtype(itype, is_itype=True)
        else:
            itype_str = itype
        validated_kwargs = {
            "name": name,
            "description": description,
            "type": type,
            "is_type": is_type,
            "dtype": dtype,
            "otype": otype,
            "n": n_features,
            "itype": itype_str,
            "minimal_set": minimal_set,
            "ordered_set": ordered_set,
            "maximal_set": maximal_set,
        }
        n_features_default = -1
        coerce_dtype_default = False
        if coerce_dtype:
            validated_kwargs["_aux"] = {"af": {"0": coerce_dtype}}
        if slots:
            list_for_hashing = [component.hash for component in slots.values()]
        else:
            HASH_CODE = {
                "dtype": "a",
                "itype": "b",
                "minimal_set": "c",
                "ordered_set": "d",
                "maximal_set": "e",
                "flexible": "f",
                "coerce_dtype": "g",
                "n": "h",
                "optional": "i",
                "features_hash": "j",
            }
            # we do not want pure informational annotations like otype, name, type, is_type, otype to be part of the hash
            hash_args = ["dtype", "itype", "minimal_set", "ordered_set", "maximal_set"]
            list_for_hashing = [
                f"{HASH_CODE[arg]}={validated_kwargs[arg]}"
                for arg in hash_args
                if validated_kwargs[arg] is not None
            ]
            # only include in hash if not default so that it's backward compatible with records for which flexible was never set
            if flexible != flexible_default:
                list_for_hashing.append(f"{HASH_CODE['flexible']}={flexible}")
            if coerce_dtype != coerce_dtype_default:
                list_for_hashing.append(f"{HASH_CODE['coerce_dtype']}={coerce_dtype}")
            if n_features != n_features_default:
                list_for_hashing.append(f"{HASH_CODE['n']}={n_features}")
            if features:
                if optional_features:
                    feature_list_for_hashing = [
                        feature.uid
                        if feature not in set(optional_features)
                        else f"{feature.uid}({HASH_CODE['optional']})"
                        for feature in features
                    ]
                else:
                    feature_list_for_hashing = [feature.uid for feature in features]
                # order matters if ordered_set is True
                if ordered_set:
                    features_hash = hash_string(":".join(feature_list_for_hashing))
                else:
                    features_hash = hash_string(
                        ":".join(sorted(feature_list_for_hashing))
                    )
                list_for_hashing.append(f"{HASH_CODE['features_hash']}={features_hash}")
        self._list_for_hashing = sorted(list_for_hashing)
        schema_hash = hash_string(":".join(self._list_for_hashing))
        validated_kwargs["hash"] = schema_hash
        return (
            features,
            validated_kwargs,
            optional_features,
            features_registry,
            flexible,
        )

    @classmethod
    def from_values(  # type: ignore
        cls,
        values: ListLike,
        field: FieldAttr = Feature.name,
        type: str | None = None,
        name: str | None = None,
        mute: bool = False,
        organism: SQLRecord | str | None = None,
        source: SQLRecord | None = None,
        raise_validation_error: bool = True,
    ) -> Schema:
        """Create feature set for validated features.

        Args:
            values: A list of values, like feature names or ids.
            field: The field of a reference registry to map values.
            type: The simple type.
                Defaults to `None` if reference registry is :class:`~lamindb.Feature`,
                defaults to `"float"` otherwise.
            name: A name.
            organism: An organism to resolve gene mapping.
            source: A public ontology to resolve feature identifier mapping.
            raise_validation_error: Whether to raise a validation error if some values are not valid.

        Raises:
            ValidationError: If some values are not valid.

        Example:

            ::

                import lamindb as ln
                import bionty as bt

                features = [ln.Feature(name=feat, dtype="str").save() for feat in ["feat11", "feat21"]]
                schema = ln.Schema.from_values(features)

                genes = ["ENSG00000139618", "ENSG00000198786"]
                schema = ln.Schema.from_values(features, bt.Gene.ensembl_gene_id, "float")
        """
        if not isinstance(field, FieldAttr):
            raise TypeError(
                "Argument `field` must be a SQLRecord field, e.g., `Feature.name`"
            )
        if len(values) == 0:
            raise ValueError("Provide a list of at least one value")
        if isinstance(values, DICT_KEYS_TYPE):
            values = list(values)
        registry = field.field.model
        if registry != Feature and type is None:
            type = NUMBER_TYPE
            logger.debug("setting feature set to 'number'")
        validated = registry.validate(values, field=field, mute=mute, organism=organism)
        values_array = np.array(values)
        validated_values = values_array[validated]
        if validated.sum() != len(values):
            not_validated_values = values_array[~validated]
            msg = (
                f"These values could not be validated: {not_validated_values.tolist()}\n"
                f"If there are no typos, add them to their registry: {registry.__name__}"
            )
            if raise_validation_error:
                raise ValidationError(msg)
            elif len(validated_values) == 0:
                return None  # temporarily return None here
        validated_features = registry.from_values(
            validated_values,
            field=field,
            organism=organism,
            source=source,
        )
        schema = Schema(
            features=validated_features,
            name=name,
            dtype=get_type_str(type),
        )
        return schema

    @classmethod
    def from_df(
        cls,
        df: pd.DataFrame,
        field: FieldAttr = Feature.name,
        name: str | None = None,
        mute: bool = False,
        organism: SQLRecord | str | None = None,
        source: SQLRecord | None = None,
    ) -> Schema | None:
        """Create schema for valid columns."""
        registry = field.field.model
        validated = registry.validate(
            df.columns, field=field, mute=mute, organism=organism
        )
        if validated.sum() == 0:
            if not mute:
                logger.warning("no validated features, skip creating schema")
            return None
        if registry == Feature:
            validated_features = Feature.from_values(  # type: ignore
                df.columns, field=field, organism=organism
            )
            schema = Schema(
                list(validated_features), name=name, dtype=None, otype="DataFrame"
            )
        else:
            dtypes = [col.dtype for (_, col) in df.loc[:, validated].items()]
            if len(set(dtypes)) != 1:
                raise ValueError(f"data types are heterogeneous: {set(dtypes)}")
            dtype = serialize_pandas_dtype(dtypes[0])
            validated_features = registry.from_values(
                df.columns[validated],
                field=field,
                organism=organism,
                source=source,
            )
            schema = Schema(
                features=list(validated_features),
                name=name,
                dtype=get_type_str(dtype),
            )
        return schema

    def save(self, *args, **kwargs) -> Schema:
        """Save."""
        from .save import bulk_create

        if self.pk is not None:
            features = (
                self._features[1]
                if hasattr(self, "_features")
                else (self.members.list() if self.members.exists() else [])
            )
            _, validated_kwargs, _, _, _ = self._validate_kwargs_calculate_hash(
                features=features,  # type: ignore
                index=None,  # need to pass None here as otherwise counting double
                slots=self.slots,
                name=self.name,
                description=self.description,
                itype=self.itype,
                flexible=self.flexible,
                type=self.type,
                is_type=self.is_type,
                otype=self.otype,
                dtype=self.dtype,
                minimal_set=self.minimal_set,
                ordered_set=self.ordered_set,
                maximal_set=self.maximal_set,
                coerce_dtype=self.coerce_dtype,
                n_features=self.n,
                optional_features_manual=self.optionals.get(),
            )
            if validated_kwargs["hash"] != self.hash:
                from .artifact import Artifact

                datasets = Artifact.filter(schema=self).all()
                if datasets.exists():
                    logger.warning(
                        f"you updated the schema hash and might invalidate datasets that were previously validated with this schema: {datasets.list('uid')}"
                    )
                self.hash = validated_kwargs["hash"]
                self.n = validated_kwargs["n"]
        super().save(*args, **kwargs)
        if hasattr(self, "_slots"):
            # analogous to save_schema_links in core._data.py
            # which is called to save feature sets in artifact.save()
            links = []
            for slot, component in self._slots.items():
                kwargs = {
                    "composite_id": self.id,
                    "component_id": component.id,
                    "slot": slot,
                }
                links.append(Schema.components.through(**kwargs))
            bulk_create(links, ignore_conflicts=True)
            delattr(self, "_slots")
        if hasattr(self, "_features"):
            assert self.n > 0  # noqa: S101
            using: bool | None = kwargs.pop("using", None)
            related_name, records = self._features
            # only the following method preserves the order
            # .set() does not preserve the order but orders by
            # the feature primary key
            through_model = getattr(self, related_name).through
            related_model_split = parse_cat_dtype(self.itype, is_itype=True)[
                "registry_str"
            ].split(".")
            if len(related_model_split) == 1:
                related_field = related_model_split[0].lower()
            else:
                related_field = related_model_split[1].lower()
            related_field_id = f"{related_field}_id"
            links = [
                through_model(**{"schema_id": self.id, related_field_id: record.id})
                for record in records
            ]
            through_model.objects.using(using).bulk_create(links, ignore_conflicts=True)
            delattr(self, "_features")
        return self

    @property
    def members(self) -> QuerySet:
        """A queryset for the individual records in the feature set underlying the schema.

        Unlike `schema.features`, `schema.genes`, `schema.proteins`, etc., this queryset is ordered and
        doesn't require knowledge of the entity.
        """
        if self._state.adding:
            # this should return a queryset and not a list...
            # need to fix this
            return self._features[1]
        if self.itype == "Composite":
            return Feature.objects.none()
        related_name = self._get_related_name()
        if related_name is None:
            related_name = "features"
        return self.__getattribute__(related_name).order_by("links_schema__id")

    @property
    def coerce_dtype(self) -> bool:
        """Whether dtypes should be coerced during validation.

        For example, a `objects`-dtyped pandas column can be coerced to `categorical` and would pass validation if this is true.
        """
        if self._aux is not None and "af" in self._aux and "0" in self._aux["af"]:  # type: ignore
            return self._aux["af"]["0"]  # type: ignore
        else:
            return False

    @coerce_dtype.setter
    def coerce_dtype(self, value: bool) -> None:
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["0"] = value

    @property
    def flexible(self) -> bool:
        """Indicates how to handle validation and annotation in case features are not defined.

        Examples:

            Make a rigid schema flexible::

                schema = ln.Schema.get(name="my_schema")
                schema.flexible = True
                schema.save()

            During schema creation::

                # if you're not passing features but just defining the itype, defaults to flexible = True
                schema = ln.Schema(itype=ln.Feature).save()
                assert not schema.flexible

                # if you're passing features, defaults to flexible = False
                schema = ln.Schema(
                    features=[ln.Feature(name="my_required_feature", dtype=int).save()],
                )
                assert not schema.flexible

                # you can also validate & annotate features in addition to those that you're explicitly defining:
                schema = ln.Schema(
                    features=[ln.Feature(name="my_required_feature", dtype=int).save()],
                    flexible=True,
                )
                assert schema.flexible

        """
        if self._aux is not None and "af" in self._aux and "2" in self._aux["af"]:  # type: ignore
            return self._aux["af"]["2"]  # type: ignore
        else:
            return (
                self.n < 0
            )  # is the flexible default, needed for backward compat if flexible was never set

    @flexible.setter
    def flexible(self, value: bool) -> None:
        self._aux = self._aux or {}
        self._aux.setdefault("af", {})["2"] = value

    @property
    def index(self) -> None | Feature:
        """The feature configured to act as index.

        To unset it, set `schema.index` to `None`.
        """
        if self._index_feature_uid is None:
            return None
        else:
            return self.features.get(uid=self._index_feature_uid)

    @index.setter
    def index(self, value: None | Feature) -> None:
        if value is None:
            current_index = self.index
            self.features.remove(current_index)
            self._index_feature_uid = value
        else:
            self.features.add(value)
            self._index_feature_uid = value.uid

    @property
    def _index_feature_uid(self) -> None | str:
        """The uid of the index feature."""
        if self._aux is not None and "af" in self._aux and "3" in self._aux["af"]:
            return self._aux["af"]["3"]
        else:
            return None

    @_index_feature_uid.setter
    def _index_feature_uid(self, value: str | None) -> None:
        self._aux = self._aux or {}
        if value is None:
            self._aux.get("af", {}).pop("3")
        else:
            self._aux.setdefault("af", {})["3"] = value

    @property
    def slots(self) -> dict[str, Schema]:
        """Slots.

        Examples:

            ::

                # define composite schema
                anndata_schema = ln.Schema(
                    name="small_dataset1_anndata_schema",
                    otype="AnnData",
                    slots={"obs": obs_schema, "var": var_schema},
                ).save()

                # access slots
                anndata_schema.slots
                # {'obs': <Schema: obs_schema>, 'var': <Schema: var_schema>}
        """
        if hasattr(self, "_slots"):
            return self._slots
        if self.itype == "Composite":
            self._slots = {
                link.slot: link.component
                for link in self.components.through.filter(composite_id=self.id).all()
            }
            return self._slots
        return {}

    @property
    def optionals(self) -> SchemaOptionals:
        """Manage optional features.

        Example:

            ::

                # a schema with optional "sample_name"
                schema_optional_sample_name = ln.Schema(
                    features=[
                        ln.Feature(name="sample_id", dtype=str).save(),  # required
                        ln.Feature(name="sample_name", dtype=str).save().with_config(optional=True),  # optional
                    ],
                ).save()

                # raise ValidationError since `sample_id` is required
                ln.curators.DataFrameCurator(
                    pd.DataFrame(
                        {
                        "sample_name": ["Sample 1", "Sample 2"],
                        }
                    ),
                    schema=schema_optional_sample_name).validate()
                )

                # passes because an optional column is missing
                ln.curators.DataFrameCurator(
                    pd.DataFrame(
                        {
                        "sample_id": ["sample1", "sample2"],
                        }
                    ),
                    schema=schema_optional_sample_name).validate()
                )
        """
        return SchemaOptionals(self)

    def describe(self, return_str=False) -> None | str:
        """Describe schema."""
        message = str(self)
        # display slots for composite schemas
        if self.itype == "Composite":
            message + "\nslots:"
            for slot, schema in self.slots.items():
                message += f"\n    {slot}: " + str(schema)
        else:
            tree = describe_schema(self)
            return format_rich_tree(
                tree, fallback="no linked features", return_str=return_str
            )
        if return_str:
            return message
        else:
            print(message)
            return None


def get_type_str(dtype: str | None) -> str | None:
    if dtype is not None:
        type_str = dtype.__name__ if not isinstance(dtype, str) else dtype  # type: ignore
    else:
        type_str = None
    return type_str


def _get_related_name(self: Schema) -> str:
    related_models = dict_related_model_to_related_name(self, instance=self._state.db)
    related_name = related_models.get(
        parse_cat_dtype(self.itype, is_itype=True)["registry_str"]
    )
    return related_name


class SchemaFeature(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    schema: Schema = ForeignKey(Schema, CASCADE, related_name="links_feature")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_schema")

    class Meta:
        unique_together = ("schema", "feature")


class ArtifactSchema(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey("Artifact", CASCADE, related_name="_links_schema")
    schema: Schema = ForeignKey(Schema, PROTECT, related_name="_links_artifact")
    slot: str | None = CharField(null=True)
    feature_ref_is_semantic: bool | None = BooleanField(null=True)

    class Meta:
        unique_together = (("artifact", "schema"), ("artifact", "slot"))


class SchemaComponent(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    composite: Schema = ForeignKey(Schema, CASCADE, related_name="links_composite")
    component: Schema = ForeignKey(Schema, PROTECT, related_name="links_component")
    slot: str | None = CharField(null=True)

    class Meta:
        unique_together = (("composite", "slot", "component"), ("composite", "slot"))


Schema._get_related_name = _get_related_name
# excluded on docs via
# https://github.com/laminlabs/lndocs/blob/8c1963de65445107ea69b3fd59354c3828e067d1/lndocs/lamin_sphinx/__init__.py#L584-L588
delattr(Schema, "validated_by")  # we don't want to expose these
delattr(Schema, "validated_by_id")  # we don't want to expose these
delattr(Schema, "validated_schemas")  # we don't want to expose these
delattr(Schema, "composite")  # we don't want to expose these
delattr(Schema, "composite_id")  # we don't want to expose these
