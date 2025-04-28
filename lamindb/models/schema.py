from __future__ import annotations

from typing import TYPE_CHECKING, Any, Type, overload

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
from lamindb.errors import FieldValidationError, InvalidArgument
from lamindb.models.feature import parse_cat_dtype

from ..errors import ValidationError
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
from .record import (
    BasicRecord,
    LinkORM,
    Record,
    Registry,
    _get_record_kwargs,
    init_self_from_db,
    update_attributes,
)
from .run import Param, TracksRun, TracksUpdates

if TYPE_CHECKING:
    from collections.abc import Iterable

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


def get_features_config(
    features: list[Record] | tuple[Record, dict],
) -> tuple[list[Record], list[tuple[Record, dict]]]:
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
        """Set the optional features."""
        if not isinstance(features, list) or not all(
            isinstance(f, Feature) for f in features
        ):
            raise TypeError("features must be a list of Feature records!")
        self.schema._aux = self.schema._aux or {}
        if len(features) > 0:
            self.schema._aux.setdefault("af", {})["1"] = [f.uid for f in features]

    def add(self, features: Feature | list[Feature]) -> None:
        """Add feature to the optional features."""
        self.schema._aux = self.schema._aux or {}
        if not isinstance(features, list):
            features = [features]
        if not all(isinstance(f, Feature) for f in features):
            raise TypeError("features must be a list of Feature records!")
        if len(features) > 0:
            if "1" not in self.schema._aux.setdefault("af", {}):
                self.set(features)
            self.schema._aux.setdefault("af", {})["1"].extend([f.uid for f in features])


class Schema(Record, CanCurate, TracksRun):
    """Schemas.

    A simple schema is a feature set such as the set of columns of a `DataFrame`.

    A composite schema has multiple components, e.g., for an `AnnData`, one schema for `obs` and another one for `var`.

    A schema can also merely define abstract constraints or instructions for dataset validation & annotation.

    Args:
        features: `Iterable[Record] | None = None` An iterable of :class:`~lamindb.Feature`
            records to hash, e.g., `[Feature(...), Feature(...)]`. Is turned into
            a set upon instantiation. If you'd like to pass values, use
            :meth:`~lamindb.Schema.from_values` or
            :meth:`~lamindb.Schema.from_df`.
        index: `Feature | None = None` A :class:`~lamindb.Feature` record to validate an index of a `DataFrame`.
        components: `dict[str, Schema] | None = None` A dictionary mapping slot names to
            components. A component is itself a :class:`~lamindb.Schema` object.
        name: `str | None = None` A name.
        description: `str | None = None` A description.
        itype: `str | None = None` The feature identifier type (e.g. :class:`~lamindb.Feature`, :class:`~bionty.Gene`, ...).
        flexible: `bool | None = None` Whether to include any feature of the same `itype` in validation
            and annotation. If no features are passed, defaults to `True`, otherwise to `False`.
        type: `Schema | None = None` A type.
        is_type: `bool = False` Distinguish types from instances of the type.
        otype: `str | None = None` An object type to define the structure of a composite schema.
        dtype: `str | None = None` The simple type. Defaults to
            `None` for sets of :class:`~lamindb.Feature` records.
            Otherwise defaults to `"num"` (e.g., for sets of :class:`~bionty.Gene`).
        minimal_set: `bool = True` Whether all passed features are to be considered required by default.
            See :attr:`~lamindb.Schema.optionals` for more-fine-grained control.
        ordered_set: `bool = False` Whether features are required to be ordered.
        maximal_set: `bool = False` If `True`, no additional features are allowed.
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

            # a schema that requires a single feature and accepts any other features with valid feature names
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

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
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
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=20)
    """A universal id (hash of the set of feature values)."""
    name: str | None = CharField(max_length=150, null=True, db_index=True)
    """A name."""
    description: str | None = CharField(null=True, db_index=True)
    """A description."""
    n: int = IntegerField()
    """Number of features in the schema."""
    itype: str | None = CharField(
        max_length=120, db_index=True, null=True, editable=False
    )
    """A registry that stores feature identifiers used in this schema, e.g., `'Feature'` or `'bionty.Gene'`.

    Depending on `itype`, `.members` stores, e.g., `Feature` or `bionty.Gene` records.
    """
    type: Schema | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of schema.

    Allows to group schemas by type, e.g., all meassurements evaluating gene expression vs. protein expression vs. multi modal.

    You can define types via `ln.Schema(name="ProteinPanel", is_type=True)`.

    Here are a few more examples for type names: `'ExpressionPanel'`, `'ProteinPanel'`, `'Multimodal'`, `'Metadata'`, `'Embedding'`.
    """
    records: Schema
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
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
    params: Param
    """The params contained in the schema."""
    artifacts: Artifact
    """The artifacts that measure a feature set that matches this schema."""
    validated_artifacts: Artifact
    """The artifacts that were validated against this schema with a :class:`~lamindb.curators.Curator`."""
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
    # """The schemas that were validated against this schema with a :class:`~lamindb.curators.Curator`."""
    composite: Schema | None = ForeignKey(
        "self", PROTECT, related_name="+", default=None, null=True
    )
    # The legacy foreign key
    slot: str | None = CharField(max_length=100, db_index=True, null=True)
    # The legacy slot

    @overload
    def __init__(
        self,
        features: Iterable[Record] | None = None,
        index: Feature | None = None,
        components: dict[str, Schema] | None = None,
        name: str | None = None,
        description: str | None = None,
        dtype: str | Type[int | float | str] | None = None,  # noqa
        itype: str | Registry | FieldAttr | None = None,
        type: Schema | None = None,
        is_type: bool = False,
        otype: str | None = None,
        ordered_set: bool = False,
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

        features: Iterable[Record] | None = (
            args[0] if args else kwargs.pop("features", [])
        )
        index: Feature | None = kwargs.pop("index", None)
        components: dict[str, Schema] = kwargs.pop("components", {})
        name: str | None = kwargs.pop("name", None)
        description: str | None = kwargs.pop("description", None)
        itype: str | Record | DeferredAttribute | None = kwargs.pop("itype", None)
        flexible: bool | None = kwargs.pop("flexible", None)
        type: Feature | None = kwargs.pop("type", None)
        is_type: bool = kwargs.pop("is_type", False)
        otype: str | None = kwargs.pop("otype", None)
        dtype: str | None = kwargs.pop("dtype", None)
        minimal_set: bool = kwargs.pop("minimal_set", True)
        ordered_set: bool = kwargs.pop("ordered_set", False)
        maximal_set: bool = kwargs.pop("maximal_set", False)
        coerce_dtype: bool | None = kwargs.pop("coerce_dtype", None)
        n_features: int | None = kwargs.pop("n", None)
        optional_features = []

        if kwargs:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Schema)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        optional_features = []
        if itype is not None:
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
            assert n_features is None, "do not pass `n` if features are passed"  # noqa: S101
            n_features = len(features)
            if features_registry == Feature:
                optional_features = [
                    config[0] for config in configs if config[1].get("optional")
                ]
        elif n_features is None:
            n_features = -1
        if dtype is None:
            dtype = None if itype is not None and itype == "Feature" else NUMBER_TYPE
        else:
            dtype = get_type_str(dtype)
        if flexible is None:
            flexible = n_features < 0
        components: dict[str, Schema]
        if components:
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
        if coerce_dtype:
            validated_kwargs["_aux"] = {"af": {"0": coerce_dtype}}
        if components:
            hash = hash_set({component.hash for component in components.values()})
        else:
            # we do not want pure informational annotations like otype, name, type, is_type, otype to be part of the hash
            hash_args = ["dtype", "itype", "minimal_set", "ordered_set", "maximal_set"]
            union_set = {
                str(validated_kwargs[arg])
                for arg in hash_args
                if validated_kwargs[arg] is not None
            }
            if flexible != n_features < 0:
                union_set.add(f"flexible:{flexible}")
            if features:
                union_set = union_set.union({feature.uid for feature in features})
            if optional_features:
                union_set = union_set.union(
                    {f"optional:{feature.uid}" for feature in optional_features}
                )
            hash = hash_set(union_set)
        validated_kwargs["hash"] = hash
        schema = Schema.filter(hash=hash).one_or_none()
        if schema is not None:
            logger.important(f"returning existing schema with same hash: {schema}")
            init_self_from_db(self, schema)
            update_attributes(self, validated_kwargs)
            self.optionals.set(optional_features)
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
            self._slots = components
        validated_kwargs["uid"] = ids.base62_20()
        super().__init__(**validated_kwargs)
        self.optionals.set(optional_features)
        self.flexible = flexible
        if index is not None:
            self._index_feature_uid = index.uid

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
        df: pd.DataFrame,
        field: FieldAttr = Feature.name,
        name: str | None = None,
        mute: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
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
                validated_features, name=name, dtype=None, otype="DataFrame"
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
                features=validated_features,
                name=name,
                dtype=get_type_str(dtype),
            )
        return schema

    def save(self, *args, **kwargs) -> Schema:
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
    def members(self) -> QuerySet:
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
            return self.n < 0

    @flexible.setter
    def flexible(self, value: bool) -> None:
        if value != (self.n < 0):
            self._aux = self._aux or {}
            self._aux.setdefault("af", {})["2"] = value

    @property
    def index(self) -> None | Feature:
        """The feature configured to act as index."""
        if self._index_feature_uid is None:
            return None
        else:
            return self.features.get(uid=self._index_feature_uid)

    @property
    def _index_feature_uid(self) -> None | str:
        """The uid of the index feature."""
        if self._aux is not None and "af" in self._aux and "3" in self._aux["af"]:
            return self._aux["af"]["3"]
        else:
            return None

    @_index_feature_uid.setter
    def _index_feature_uid(self, value: str) -> None:
        self._aux = self._aux or {}
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
                    components={"obs": obs_schema, "var": var_schema},
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
    related_name = related_models.get(parse_cat_dtype(self.itype)["registry_str"])
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
    artifact: Artifact = ForeignKey("Artifact", CASCADE, related_name="_links_schema")
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
