from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Optional, overload

import numpy as np
from django.db import models
from django.db.models import CASCADE, PROTECT, ManyToManyField
from lamin_utils import logger
from lamindb_setup.core.hashing import HASH_LENGTH, hash_set

from lamindb.base import ids
from lamindb.base.fields import (
    BooleanField,
    CharField,
    ForeignKey,
    IntegerField,
    JSONField,
)
from lamindb.base.types import FieldAttr, ListLike
from lamindb.errors import InvalidArgument

from ..base import deprecated
from ..errors import ValidationError
from ._relations import (
    dict_related_model_to_related_name,
    get_related_name,
)
from .base import LinkORM, TracksRun, TracksUpdates
from .can_curate import CanCurate
from .core import Param
from .feature import (
    Feature,
    convert_pandas_dtype_to_lamin_dtype,
    get_dtype_str_from_dtype,
)
from .record import (
    BasicRecord,
    Record,
    Registry,
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


def validate_features(features: list[Record]) -> Record:
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
    if not isinstance(features[0], Record):
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


class Schema(Record, CanCurate, TracksRun):
    """Schemas / feature sets.

    A simple schema is just a set of columns in a `DataFrame`, a "feature set".

    A composite schema has multiple components, e.g. for an `AnnData`, each a feature set for `obs` and `var`.

    Args:
        features: `Iterable[Record] | None = None` An iterable of :class:`~lamindb.Feature`
            records to hash, e.g., `[Feature(...), Feature(...)]`. Is turned into
            a set upon instantiation. If you'd like to pass values, use
            :meth:`~lamindb.Schema.from_values` or
            :meth:`~lamindb.Schema.from_df`.
        components: `dict[str, Schema] | None = None` A dictionary mapping component names to
            their corresponding :class:`~lamindb.Schema` objects for composite schemas.
        name: `str | None = None` A name.
        description: `str | None = None` A description.
        dtype: `str | None = None` The simple type. Defaults to
            `None` for sets of :class:`~lamindb.Feature` records.
            Otherwise defaults to `"num"` (e.g., for sets of :class:`~bionty.Gene`).
        itype: `str | None = None` The feature identifier type (e.g. :class:`~lamindb.Feature`, :class:`~bionty.Gene`, ...).
        type: `Schema | None = None` A type.
        is_type: `bool = False` Distinguish types from instances of the type.
        otype: `str | None = None` An object type to define the structure of a composite schema.
        minimal_set: `bool = True` Whether the schema contains a minimal set of linked features.
        ordered_set: `bool = False` Whether features are required to be ordered.
        maximal_set: `bool = False` If `True`, no additional features are allowed.
        slot: `str | None = None` The slot name when this schema is used as a component in a
            composite schema.
        coerce_dtype: `bool = False` When True, attempts to coerce values to the specified dtype
            during validation, see :attr:`~lamindb.Schema.coerce_dtype`.

    .. dropdown:: Why does LaminDB model schemas, not just features?

        1. Performance: Imagine you measure the same panel of 20k transcripts in
           1M samples. By modeling the panel as a feature set, you can link all
           your artifacts against one feature set and only need to store 1M
           instead of 1M x 20k = 20B links.
        2. Interpretation: Model protein panels, gene panels, etc.
        3. Data integration: Feature sets provide the information that determines whether two datasets can be meaningfully concatenated.

        These reasons do not hold for label sets. Hence, LaminDB does not model label sets.

    Note:

        A feature set can be identified by the `hash` of its feature uids.
        It's stored in the `.hash` field.

        A `slot` provides a string key to access feature sets. For instance, for the schema of an
        `AnnData` object, it would be `'obs'` for `adata.obs`.

    See Also:
        :meth:`~lamindb.Schema.from_values`
            Create from values.
        :meth:`~lamindb.Schema.from_df`
            Create from dataframe columns.

    Examples:

        Create a schema (feature set) from df with types:

        >>> df = pd.DataFrame({"feat1": [1, 2], "feat2": [3.1, 4.2], "feat3": ["cond1", "cond2"]})
        >>> schema = ln.Schema.from_df(df)

        Create a schema (feature set) from features:

        >>> features = [ln.Feature(name=feat, dtype="float").save() for feat in ["feat1", "feat2"]]
        >>> schema = ln.Schema(features)

        Create a schema (feature set) from identifier values:

        >>> import bionty as bt
        >>> schema = ln.Schema.from_values(adata.var["ensemble_id"], Gene.ensembl_gene_id, organism="mouse").save()

    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"
    _aux_fields: dict[str, tuple[str, type]] = {
        "0": ("coerce_dtype", bool),
        "1": ("_index_feature_uid", str),
    }

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=20)
    """A universal id (hash of the set of feature values)."""
    name: str | None = CharField(max_length=150, null=True, db_index=True)
    """A name."""
    description: str | None = CharField(null=True, db_index=True)
    """A description."""
    n = IntegerField()
    """Number of features in the set."""
    dtype: str | None = CharField(max_length=64, null=True, editable=False)
    """Data type, e.g., "num", "float", "int". Is `None` for :class:`~lamindb.Feature`.

    For :class:`~lamindb.Feature`, types are expected to be heterogeneous and defined on a per-feature level.
    """
    itype: str | None = CharField(
        max_length=120, db_index=True, null=True, editable=False
    )
    """A registry that stores feature identifiers used in this schema, e.g., `'Feature'` or `'bionty.Gene'`.

    Depending on the registry, `.members` stores, e.g., `Feature` or `bionty.Gene` records.

    .. versionchanged:: 1.0.0
        Was called `registry` before.
    """
    type: Optional["Schema"] = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of schema.

    Allows to group schemas by type, e.g., all meassurements evaluating gene expression vs. protein expression vs. multi modal.

    You can define types via `ln.Schema(name="ProteinPanel", is_type=True)`.

    Here are a few more examples for type names: `'ExpressionPanel'`, `'ProteinPanel'`, `'Multimodal'`, `'Metadata'`, `'Embedding'`.
    """
    records: "Schema"
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    otype: str | None = CharField(max_length=64, db_index=True, null=True)
    """Default Python object type, e.g., DataFrame, AnnData."""
    hash: str | None = CharField(
        max_length=HASH_LENGTH, db_index=True, null=True, editable=False
    )
    """A hash of the set of feature identifiers.

    For a composite schema, the hash of hashes.
    """
    minimal_set: bool = BooleanField(default=True, db_index=True, editable=False)
    """Whether the schema contains a minimal set of linked features (default `True`).

    If `False`, no features are linked to this schema.

    If `True`, features are linked and considered as a minimally required set in validation.
    """
    ordered_set: bool = BooleanField(default=False, db_index=True, editable=False)
    """Whether features are required to be ordered (default `False`)."""
    maximal_set: bool = BooleanField(default=False, db_index=True, editable=False)
    """If `False`, additional features are allowed (default `False`).

    If `True`, the the minimal set is a maximal set and no additional features are allowed.
    """
    components: "Schema" = ManyToManyField(
        "self", through="SchemaComponent", symmetrical=False, related_name="composites"
    )
    """Components of this schema."""
    composites: "Schema"
    """The composite schemas that contains this schema as a component.

    For example, an `AnnData` composes multiple schemas: `var[DataFrameT]`, `obs[DataFrame]`, `obsm[Array]`, `uns[dict]`, etc.
    """
    features: Feature
    """The features contained in the schema."""
    params: Param
    """The params contained in the schema."""
    artifacts: "Artifact"
    """The artifacts that measure a feature set that matches this schema."""
    validated_artifacts: "Artifact"
    """The artifacts that were validated against this schema with a :class:`~lamindb.curators.Curator`."""
    projects: "Project"
    """Associated projects."""
    _curation: dict[str, Any] = JSONField(default=None, db_default=None, null=True)
    # lamindb v2
    # _itype: ContentType = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    # ""Index of the registry that stores the feature identifiers, e.g., `Feature` or `Gene`."""
    # -- the following two fields are dynamically removed from the API for now
    validated_by: Optional["Schema"] = ForeignKey(
        "self", PROTECT, related_name="validated_schemas", default=None, null=True
    )
    # """The schema that validated this schema during curation.

    # When performing validation, the schema that enforced validation is often less concrete than what is validated.

    # For instance, the set of measured features might be a superset of the minimally required set of features.
    # """
    # validated_schemas: Schema
    # """The schemas that were validated against this schema with a :class:`~lamindb.curators.Curator`."""
    composite: Optional["Schema"] = ForeignKey(
        "self", PROTECT, related_name="+", default=None, null=True
    )
    # The legacy foreign key
    slot: str | None = CharField(max_length=100, db_index=True, null=True)
    # The legacy slot

    @overload
    def __init__(
        self,
        features: Iterable[Record] | None = None,
        components: dict[str, "Schema"] | None = None,
        name: str | None = None,
        description: str | None = None,
        dtype: str | None = None,
        itype: str | Registry | FieldAttr | None = None,
        type: Optional["Schema"] = None,
        is_type: bool = False,
        otype: str | None = None,
        minimal_set: bool = True,
        ordered_set: bool = False,
        maximal_set: bool = False,
        slot: str | None = None,
        coerce_dtype: bool = False,
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

        features: Iterable[Record] | None = (
            args[0] if args else kwargs.pop("features", [])
        )
        # typing here anticipates transitioning to a ManyToMany
        # between composites and components similar to feature_sets
        # in lamindb v2
        components: dict[str, Schema] = kwargs.pop("components", {})
        name: str | None = kwargs.pop("name", None)
        description: str | None = kwargs.pop("description", None)
        dtype: str | None = kwargs.pop("dtype", None)
        itype: str | Record | DeferredAttribute | None = kwargs.pop("itype", None)
        type: Feature | None = kwargs.pop("type", None)
        is_type: bool = kwargs.pop("is_type", False)
        otype: str | None = kwargs.pop("otype", None)
        minimal_set: bool = kwargs.pop("minimal_set", True)
        ordered_set: bool = kwargs.pop("ordered_set", False)
        maximal_set: bool = kwargs.pop("maximal_set", False)
        slot: str | None = kwargs.pop("slot", None)
        coerce_dtype: bool | None = kwargs.pop("coerce_dtype", None)

        if kwargs:
            raise ValueError(
                f"Unexpected keyword arguments: {', '.join(kwargs.keys())}\n"
                "Valid arguments are: features, description, dtype, itype, type, "
                "is_type, otype, minimal_set, ordered_set, maximal_set, "
                "slot, validated_by, coerce_dtype"
            )

        if features:
            features_registry = validate_features(features)
            itype_compare = features_registry.__get_name_with_module__()
            if itype is not None:
                assert itype == itype_compare, str(itype_compare)  # noqa: S101
            else:
                itype = itype_compare
            n_features = len(features)
        else:
            n_features = -1
        if dtype is None:
            dtype = None if itype is not None and itype == "Feature" else NUMBER_TYPE
        else:
            dtype = get_type_str(dtype)
        components: dict[str, Schema]
        if components:
            itype = "Composite"
            if otype is None:
                raise InvalidArgument("Please pass otype != None for composite schemas")
        if itype is not None and not isinstance(itype, str):
            itype_str = get_dtype_str_from_dtype(itype, is_itype=True)
        else:
            itype_str = itype
        validated_kwargs = {
            "name": name,
            "description": description,
            "type": type,
            "dtype": dtype,
            "is_type": is_type,
            "otype": otype,
            "n": n_features,
            "itype": itype_str,
            "minimal_set": minimal_set,
            "ordered_set": ordered_set,
            "maximal_set": maximal_set,
        }
        if coerce_dtype:
            validated_kwargs["_aux"] = {"af": {"0": coerce_dtype}}
        if features:
            hash = hash_set({feature.uid for feature in features})
        elif components:
            hash = hash_set({component.hash for component in components.values()})
        else:
            hash = hash_set({str(value) for value in validated_kwargs.values()})
        validated_kwargs["hash"] = hash
        validated_kwargs["slot"] = slot
        schema = Schema.filter(hash=hash).one_or_none()
        if schema is not None:
            logger.important(f"returning existing schema with same hash: {schema}")
            init_self_from_db(self, schema)
            update_attributes(self, validated_kwargs)
            return None
        self._components: dict[str, Schema] = {}
        if features:
            self._features = (get_related_name(features_registry), features)  # type: ignore
        elif components:
            for slot, component in components.items():
                if component._state.adding:
                    raise InvalidArgument(
                        f"component {slot} {component} must be saved before use"
                    )
            self._components = components
        validated_kwargs["uid"] = ids.base62_20()
        super().__init__(**validated_kwargs)

    @classmethod
    def from_values(  # type: ignore
        cls,
        values: ListLike,
        field: FieldAttr = Feature.name,
        type: str | None = None,
        name: str | None = None,
        mute: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
        raise_validation_error: bool = True,
    ) -> "Schema":
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

        Examples:

            >>> features = [ln.Feature(name=feat, dtype="str").save() for feat in ["feat11", "feat21"]]
            >>> schema = ln.Schema.from_values(features)

            >>> genes = ["ENSG00000139618", "ENSG00000198786"]
            >>> schema = ln.Schema.from_values(features, bt.Gene.ensembl_gene_id, "float")
        """
        if not isinstance(field, FieldAttr):
            raise TypeError(
                "Argument `field` must be a Record field, e.g., `Feature.name`"
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
        df: "pd.DataFrame",
        field: FieldAttr = Feature.name,
        name: str | None = None,
        mute: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
    ) -> Optional["Schema"]:
        """Create feature set for validated features."""
        registry = field.field.model
        validated = registry.validate(
            df.columns, field=field, mute=mute, organism=organism
        )
        if validated.sum() == 0:
            if mute is True:
                logger.warning("no validated features, skip creating feature set")
            return None
        if registry == Feature:
            validated_features = Feature.from_values(  # type: ignore
                df.columns, field=field, organism=organism
            )
            schema = Schema(
                validated_features, name=name, dtype=None, otype="DataFrame"
            )
        else:
            dtypes = [col.dtype for (_, col) in df.loc[:, validated].items()]
            if len(set(dtypes)) != 1:
                raise ValueError(f"data types are heterogeneous: {set(dtypes)}")
            dtype = convert_pandas_dtype_to_lamin_dtype(dtypes[0])
            validated_features = registry.from_values(
                df.columns[validated],
                field=field,
                organism=organism,
                source=source,
            )
            schema = Schema(
                features=validated_features,
                name=name,
                dtype=get_type_str(dtype),
                otype="DataFrame",
            )
        return schema

    def save(self, *args, **kwargs) -> "Schema":
        """Save."""
        from .save import bulk_create

        super().save(*args, **kwargs)
        if hasattr(self, "_components"):
            # analogous to save_schema_links in core._data.py
            # which is called to save feature sets in artifact.save()
            links = []
            for slot, component in self._components.items():
                kwargs = {
                    "composite_id": self.id,
                    "component_id": component.id,
                    "slot": slot,
                }
                links.append(Schema.components.through(**kwargs))
            bulk_create(links, ignore_conflicts=True)
        if hasattr(self, "_features"):
            assert self.n > 0  # noqa: S101
            related_name, records = self._features
            # only the following method preserves the order
            # .set() does not preserve the order but orders by
            # the feature primary key
            through_model = getattr(self, related_name).through
            related_model_split = self.itype.split(".")
            if len(related_model_split) == 1:
                related_field = related_model_split[0].lower()
            else:
                related_field = related_model_split[1].lower()
            related_field_id = f"{related_field}_id"
            links = [
                through_model(**{"schema_id": self.id, related_field_id: record.id})
                for record in records
            ]
            through_model.objects.bulk_create(links, ignore_conflicts=True)
        return self

    @property
    def members(self) -> "QuerySet":
        """A queryset for the individual records of the set."""
        if self._state.adding:
            # this should return a queryset and not a list...
            # need to fix this
            return self._features[1]
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
        if self._aux is None:  # type: ignore
            self._aux = {}  # type: ignore
        if "af" not in self._aux:
            self._aux["af"] = {}
        self._aux["af"]["0"] = value

    @coerce_dtype.setter
    def coerce_dtype(self, value: bool) -> None:
        if self._aux is None:
            self._aux = {}
        if "af" not in self._aux:
            self._aux["af"] = {}
        self._aux["af"]["0"] = value

    # @property
    # def index_feature(self) -> None | Feature:
    #     # index_feature: `Record | None = None` A :class:`~lamindb.Feature` to validate the index of a `DataFrame`.
    #     """The uid of the index feature, if `index_feature` was set."""
    #     if self._index_feature_uid is None:
    #         return None
    #     else:
    #         return self.features.get(uid=self._index_feature_uid)

    # @property
    # def _index_feature_uid(self) -> None | str:
    #     """The uid of the index feature, if `index_feature` was set."""
    #     if self._aux is not None and "af" in self._aux and "1" in self._aux["af"]:
    #         return self._aux["af"]["1"]
    #     else:
    #         return None

    # @_index_feature_uid.setter
    # def _index_feature_uid(self, value: str) -> None:
    #     if self._aux is None:
    #         self._aux = {}
    #     if "af" not in self._aux:
    #         self._aux["af"] = {}
    #     self._aux["af"]["1"] = value

    @property
    @deprecated("itype")
    def registry(self) -> str:
        return self.itype

    @registry.setter
    def registry(self, value) -> None:
        self.itype = value

    def describe(self, return_str=False) -> None | str:
        """Describe schema."""
        message = str(self) + "\ncomponents:"
        for component in self.components.all():
            message += "\n    " + str(component)
        if return_str:
            return message
        else:
            print(message)
            return None

    def _get_component(self, slot: str) -> "Schema":
        return self.components.get(links_component__slot=slot)


def get_type_str(dtype: str | None) -> str | None:
    if dtype is not None:
        type_str = dtype.__name__ if not isinstance(dtype, str) else dtype  # type: ignore
    else:
        type_str = None
    return type_str


def _get_related_name(self: Schema) -> str:
    related_models = dict_related_model_to_related_name(self, instance=self._state.db)
    related_name = related_models.get(self.itype)
    return related_name


class SchemaFeature(BasicRecord, LinkORM):
    id: int = models.BigAutoField(primary_key=True)
    schema: Schema = ForeignKey(Schema, CASCADE, related_name="links_feature")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_schema")

    class Meta:
        unique_together = ("schema", "feature")


class SchemaParam(BasicRecord, LinkORM):
    id: int = models.BigAutoField(primary_key=True)
    schema: Schema = ForeignKey(Schema, CASCADE, related_name="+")
    param: Param = ForeignKey(Param, PROTECT, related_name="+")

    class Meta:
        unique_together = ("schema", "param")


class ArtifactSchema(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: "Artifact" = ForeignKey("Artifact", CASCADE, related_name="_links_schema")
    schema: Schema = ForeignKey(Schema, PROTECT, related_name="_links_artifact")
    slot: str | None = CharField(null=True)
    feature_ref_is_semantic: bool | None = BooleanField(null=True)

    class Meta:
        unique_together = (("artifact", "schema"), ("artifact", "slot"))


class SchemaComponent(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    composite: Schema = ForeignKey(Schema, CASCADE, related_name="links_composite")
    component: Schema = ForeignKey(Schema, PROTECT, related_name="links_component")
    slot: str | None = CharField(null=True)

    class Meta:
        unique_together = (("composite", "component"), ("composite", "slot"))


Schema._get_related_name = _get_related_name
# excluded on docs via
# https://github.com/laminlabs/lndocs/blob/8c1963de65445107ea69b3fd59354c3828e067d1/lndocs/lamin_sphinx/__init__.py#L584-L588
delattr(Schema, "validated_by")  # we don't want to expose these
delattr(Schema, "validated_by_id")  # we don't want to expose these
delattr(Schema, "validated_schemas")  # we don't want to expose these
delattr(Schema, "composite")  # we don't want to expose these
delattr(Schema, "composite_id")  # we don't want to expose these
