from __future__ import annotations

import builtins
import inspect
import re
import sys
from collections import defaultdict
from functools import reduce
from itertools import chain
from pathlib import PurePosixPath
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NamedTuple,
    TypeVar,
    Union,
    overload,
)

import dj_database_url
import lamindb_setup as ln_setup
from django.core.exceptions import ValidationError as DjangoValidationError
from django.db import IntegrityError, connections, models, transaction
from django.db.models import (
    CASCADE,
    PROTECT,
    Field,
    IntegerField,
    Manager,
    Q,
    QuerySet,
    Value,
)
from django.db.models.base import ModelBase
from django.db.models.fields.related import (
    ManyToManyField,
    ManyToManyRel,
    ManyToOneRel,
)
from django.db.models.functions import Cast, Coalesce
from django.db.models.lookups import (
    Contains,
    Exact,
    IContains,
    IExact,
    IRegex,
    IStartsWith,
    Regex,
    StartsWith,
)
from lamin_utils import colors, logger
from lamin_utils._lookup import Lookup
from lamindb_setup import settings as setup_settings
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
    update_db_using_local,
)
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance_hub
from lamindb_setup.core._settings_store import instance_settings_file
from lamindb_setup.core.upath import extract_suffix_from_path

from lamindb.base import deprecated
from lamindb.base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    TextField,
)
from lamindb.base.types import FieldAttr, StrField
from lamindb.errors import FieldValidationError

from ..errors import (
    InvalidArgument,
    RecordNameChangeIntegrityError,
    ValidationError,
)
from ._is_versioned import IsVersioned

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd

    from .artifact import Artifact
    from .run import Run, User
    from .transform import Transform


T = TypeVar("T", bound="Record")
IPYTHON = getattr(builtins, "__IPYTHON__", False)


# -------------------------------------------------------------------------------------
# A note on required fields at the Record level
#
# As Django does most of its validation on the Form-level, it doesn't offer functionality
# for validating the integrity of an Record object upon instantation (similar to pydantic)
#
# For required fields, we define them as commonly done on the SQL level together
# with a validator in Record (validate_required_fields)
#
# This goes against the Django convention, but goes with the SQLModel convention
# (Optional fields can be null on the SQL level, non-optional fields cannot)
#
# Due to Django's convention where CharFieldAttr has pre-configured (null=False, default=""), marking
# a required field necessitates passing `default=None`. Without the validator it would trigger
# an error at the SQL-level, with it, it triggers it at instantiation

# -------------------------------------------------------------------------------------
# A note on class and instance methods of core Record
#
# All of these are defined and tested within lamindb, in files starting with _{orm_name}.py

# -------------------------------------------------------------------------------------
# A note on maximal lengths of char fields
#
# 100 characters:
#     "Raindrops pitter-pattered on the windowpane, blurring the"
#     "city lights outside, curled up with a mug."
# A good maximal length for a name (title).
#
# 150 characters: We choose this for name maximal length because some users like long names.
#
# 255 characters:
#     "In creating a precise 255-character paragraph, one engages in"
#     "a dance of words, where clarity meets brevity. Every syllable counts,"
#     "illustrating the skill in compact expression, ensuring the essence of the"
#     "message shines through within the exacting limit."
# This is a good maximal length for a description field.


class LinkORM:
    pass


def deferred_attribute__repr__(self):
    return f"FieldAttr({self.field.model.__name__}.{self.field.name})"


FieldAttr.__repr__ = deferred_attribute__repr__  # type: ignore


class ValidateFields:
    pass


def is_approx_pascal_case(s):
    """Check if the last component of a dotted string is in PascalCase.

    Args:
        s (str): The string to check

    Returns:
        bool: True if the last component is in PascalCase

    Raises:
        ValueError: If the last component doesn't start with a capital letter
    """
    if "[" in s:  # this is because we allow types of form 'script[test_script.py]'
        return True
    last_component = s.split(".")[-1]

    if not last_component[0].isupper():
        raise ValueError(
            f"'{last_component}' should start with a capital letter given you're defining a type"
        )

    return True


def init_self_from_db(self: Record, existing_record: Record):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def update_attributes(record: Record, attributes: dict[str, str]):
    for key, value in attributes.items():
        if (
            getattr(record, key) != value
            and value is not None
            and key != "dtype"
            and key != "_aux"
        ):
            logger.warning(f"updated {key} from {getattr(record, key)} to {value}")
            setattr(record, key, value)


def validate_literal_fields(record: Record, kwargs) -> None:
    """Validate all Literal type fields in a record.

    Args:
        record: record being validated

    Raises:
        ValidationError: If any field value is not in its Literal's allowed values
    """
    if isinstance(record, LinkORM):
        return None
    if record.__class__.__name__ in "Feature":
        return None
    from lamindb.base.types import Dtype, TransformType

    types = {
        "TransformType": TransformType,
        "ArtifactKind": Dtype,
        "Dtype": Dtype,
    }
    errors = {}
    annotations = getattr(record.__class__, "__annotations__", {})
    for field_name, annotation in annotations.items():
        if field_name not in kwargs or kwargs[field_name] is None:
            continue
        value = kwargs[field_name]
        if str(annotation) in types:
            annotation = types[annotation]
        if not hasattr(annotation, "__origin__"):
            continue
        literal_type = annotation if annotation.__origin__ is Literal else None
        if literal_type is None:
            continue
        valid_values = set(literal_type.__args__)
        if value not in valid_values:
            errors[field_name] = (
                f"{field_name}: {colors.yellow(value)} is not a valid value"
                f"\n    → Valid values are: {colors.green(', '.join(sorted(valid_values)))}"
            )
    if errors:
        message = "\n  "
        for _, error in errors.items():
            message += error + "\n  "
        raise FieldValidationError(message)


def validate_fields(record: Record, kwargs):
    from lamindb.models import (
        Artifact,
        Collection,
        Feature,
        Param,
        Run,
        Schema,
        Transform,
        ULabel,
    )

    # validate required fields
    # a "required field" is a Django field that has `null=False, default=None`
    required_fields = {
        k.name for k in record._meta.fields if not k.null and k.default is None
    }
    required_fields_not_passed = {k: None for k in required_fields if k not in kwargs}
    kwargs.update(required_fields_not_passed)
    missing_fields = [
        k for k, v in kwargs.items() if v is None and k in required_fields
    ]
    if missing_fields:
        raise FieldValidationError(f"{missing_fields} are required.")
    # ensure the exact length of the internal uid for core entities
    if "uid" in kwargs and record.__class__ in {
        Artifact,
        Collection,
        Transform,
        Run,
        ULabel,
        Feature,
        Schema,
        Param,
    }:
        uid_max_length = record.__class__._meta.get_field(
            "uid"
        ).max_length  # triggers FieldDoesNotExist
        if len(kwargs["uid"]) != uid_max_length:  # triggers KeyError
            raise ValidationError(
                f"`uid` must be exactly {uid_max_length} characters long, got {len(kwargs['uid'])}."
            )
    # validate is_type
    if "is_type" in kwargs and "name" in kwargs and kwargs["is_type"]:
        if kwargs["name"].endswith("s"):
            logger.warning(
                f"name '{kwargs['name']}' for type ends with 's', in case you're naming with plural, consider the singular for a type name"
            )
        is_approx_pascal_case(kwargs["name"])
    # validate literals
    validate_literal_fields(record, kwargs)


def suggest_records_with_similar_names(
    record: Record, name_field: str, kwargs
) -> Record | None:
    """Returns True if found exact match, otherwise False.

    Logs similar matches if found.
    """
    if kwargs.get(name_field) is None or not isinstance(kwargs.get(name_field), str):
        return None
    # need to perform an additional request to find the exact match
    # previously, this was inferred from the truncated/fuzzy search below
    # but this isn't reliable: https://laminlabs.slack.com/archives/C04FPE8V01W/p1737812808563409
    # the below needs to be .first() because there might be multiple records with the same
    # name field in case the record is versioned (e.g. for Transform key)
    exact_match = record.__class__.filter(**{name_field: kwargs[name_field]}).first()
    if exact_match is not None:
        return exact_match
    queryset = _search(
        record.__class__,
        kwargs[name_field],
        field=name_field,
        truncate_string=True,
        limit=3,
    )
    if not queryset.exists():  # empty queryset
        return None
    s, it, nots = ("", "it", "s") if len(queryset) == 1 else ("s", "one of them", "")
    msg = f"record{s} with similar {name_field}{s} exist{nots}! did you mean to load {it}?"
    if IPYTHON:
        from IPython.display import display

        from lamindb import settings

        logger.warning(f"{msg}")
        if settings._verbosity_int >= 1:
            display(queryset.df())
    else:
        logger.warning(f"{msg}\n{queryset}")
    return None


RECORD_REGISTRY_EXAMPLE = """Example::

        from lamindb import Record, fields

        # sub-classing `Record` creates a new registry
        class Experiment(Record):
            name: str = fields.CharField()

        # instantiating `Experiment` creates a record `experiment`
        experiment = Experiment(name="my experiment")

        # you can save the record to the database
        experiment.save()

        # `Experiment` refers to the registry, which you can query
        df = Experiment.filter(name__startswith="my ").df()
"""


# this is the metaclass for Record
@doc_args(RECORD_REGISTRY_EXAMPLE)
class Registry(ModelBase):
    """Metaclass for :class:`~lamindb.models.Record`.

    Each `Registry` *object* is a `Record` *class* and corresponds to a table in the metadata SQL database.

    You work with `Registry` objects whenever you use *class methods* of `Record`.

    You call any subclass of `Record` a "registry" and their objects "records". A `Record` object corresponds to a row in the SQL table.

    If you want to create a new registry, you sub-class `Record`.

    {}

    Note: `Registry` inherits from Django's `ModelBase`.
    """

    def __new__(cls, name, bases, attrs, **kwargs):
        new_class = super().__new__(cls, name, bases, attrs, **kwargs)
        return new_class

    # below creates a sensible auto-complete behavior that differs across the
    # class and instance level in Jupyter Editors it doesn't have any effect for
    # static type analyzer like pylance used in VSCode
    def __dir__(cls):
        # this is needed to bring auto-complete on the class-level back
        # https://laminlabs.slack.com/archives/C04FPE8V01W/p1717535625268849
        # Filter class attributes, excluding instance methods
        exclude_instance_methods = "sphinx" not in sys.modules
        # https://laminlabs.slack.com/archives/C04FPE8V01W/p1721134595920959

        def include_attribute(attr_name, attr_value):
            if attr_name.startswith("__"):
                return False
            if exclude_instance_methods and callable(attr_value):
                return isinstance(attr_value, (classmethod, staticmethod, type))
            return True

        # check also inherited attributes
        if hasattr(cls, "mro"):
            attrs = chain(*(c.__dict__.items() for c in cls.mro()))
        else:
            attrs = cls.__dict__.items()

        result = []
        for attr_name, attr_value in attrs:
            if attr_name not in result and include_attribute(attr_name, attr_value):
                result.append(attr_name)

        # Add non-dunder attributes from Registry
        for attr in dir(Registry):
            if not attr.startswith("__") and attr not in result:
                result.append(attr)
        return result

    def __repr__(cls) -> str:
        return registry_repr(cls)

    def lookup(
        cls,
        field: StrField | None = None,
        return_field: StrField | None = None,
    ) -> NamedTuple:
        """Return an auto-complete object for a field.

        Args:
            field: The field to look up the values for. Defaults to first string field.
            return_field: The field to return. If `None`, returns the whole record.

        Returns:
            A `NamedTuple` of lookup information of the field values with a
            dictionary converter.

        See Also:
            :meth:`~lamindb.models.Record.search`

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> bt.Gene.from_source(symbol="ADGB-DT").save()
            >>> lookup = bt.Gene.lookup()
            >>> lookup.adgb_dt
            >>> lookup_dict = lookup.dict()
            >>> lookup_dict['ADGB-DT']
            >>> lookup_by_ensembl_id = bt.Gene.lookup(field="ensembl_gene_id")
            >>> genes.ensg00000002745
            >>> lookup_return_symbols = bt.Gene.lookup(field="ensembl_gene_id", return_field="symbol")
        """
        return _lookup(cls=cls, field=field, return_field=return_field)

    def filter(cls, *queries, **expressions) -> QuerySet:
        """Query records.

        Args:
            queries: One or multiple `Q` objects.
            expressions: Fields and values passed as Django query expressions.

        Returns:
            A :class:`~lamindb.models.QuerySet`.

        See Also:
            - Guide: :doc:`docs:registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples:
            >>> ln.ULabel(name="my label").save()
            >>> ln.ULabel.filter(name__startswith="my").df()
        """
        from .query_set import QuerySet

        _using_key = None
        if "_using_key" in expressions:
            _using_key = expressions.pop("_using_key")

        return QuerySet(model=cls, using=_using_key).filter(*queries, **expressions)

    def get(
        cls: type[T],
        idlike: int | str | None = None,
        **expressions,
    ) -> T:
        """Get a single record.

        Args:
            idlike: Either a uid stub, uid or an integer id.
            expressions: Fields and values passed as Django query expressions.

        Raises:
            :exc:`docs:lamindb.errors.DoesNotExist`: In case no matching record is found.

        See Also:
            - Guide: :doc:`docs:registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples::

            ulabel = ln.ULabel.get("FvtpPJLJ")
            ulabel = ln.ULabel.get(name="my-label")
        """
        from .query_set import QuerySet

        return QuerySet(model=cls).get(idlike, **expressions)

    def df(
        cls,
        include: str | list[str] | None = None,
        features: bool | list[str] = False,
        limit: int = 100,
    ) -> pd.DataFrame:
        """Convert to `pd.DataFrame`.

        By default, shows all direct fields, except `updated_at`.

        Use arguments `include` or `feature` to include other data.

        Args:
            include: Related fields to include as columns. Takes strings of
                form `"ulabels__name"`, `"cell_types__name"`, etc. or a list
                of such strings.
            features: If `True`, map all features of the
                :class:`~lamindb.Feature` registry onto the resulting
                `DataFrame`. Only available for `Artifact`.
            limit: Maximum number of rows to display from a Pandas DataFrame.
                Defaults to 100 to reduce database load.

        Examples:

            Include the name of the creator in the `DataFrame`:

            >>> ln.ULabel.df(include="created_by__name"])

            Include display of features for `Artifact`:

            >>> df = ln.Artifact.df(features=True)
            >>> ln.view(df)  # visualize with type annotations

            Only include select features:

            >>> df = ln.Artifact.df(features=["cell_type_by_expert", "cell_type_by_model"])
        """
        query_set = cls.filter()
        if hasattr(cls, "updated_at"):
            query_set = query_set.order_by("-updated_at")
        return query_set[:limit].df(include=include, features=features)

    def search(
        cls,
        string: str,
        *,
        field: StrField | None = None,
        limit: int | None = 20,
        case_sensitive: bool = False,
    ) -> QuerySet:
        """Search.

        Args:
            string: The input string to match against the field ontology values.
            field: The field or fields to search. Search all string fields by default.
            limit: Maximum amount of top results to return.
            case_sensitive: Whether the match is case sensitive.

        Returns:
            A sorted `DataFrame` of search results with a score in column `score`.
            If `return_queryset` is `True`.  `QuerySet`.

        See Also:
            :meth:`~lamindb.models.Record.filter`
            :meth:`~lamindb.models.Record.lookup`

        Examples:
            >>> ulabels = ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")
            >>> ln.save(ulabels)
            >>> ln.ULabel.search("ULabel2")
        """
        return _search(
            cls=cls,
            string=string,
            field=field,
            limit=limit,
            case_sensitive=case_sensitive,
        )

    def using(
        cls,
        instance: str | None,
    ) -> QuerySet:
        """Use a non-default LaminDB instance.

        Args:
            instance: An instance identifier of form "account_handle/instance_name".

        Examples:
            >>> ln.ULabel.using("account_handle/instance_name").search("ULabel7", field="name")
                        uid    score
            name
            ULabel7  g7Hk9b2v  100.0
            ULabel5  t4Jm6s0q   75.0
            ULabel6  r2Xw8p1z   75.0
        """
        from .query_set import QuerySet

        # we're in the default instance
        if instance is None or instance == "default":
            return QuerySet(model=cls, using=None)
        owner, name = get_owner_name_from_identifier(instance)
        if [owner, name] == setup_settings.instance.slug.split("/"):
            return QuerySet(model=cls, using=None)

        # move on to different instances
        settings_file = instance_settings_file(name, owner)
        cache_filepath = (
            ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
        )
        if not settings_file.exists():
            result = connect_instance_hub(owner=owner, name=name)
            if isinstance(result, str):
                raise RuntimeError(
                    f"Failed to load instance {instance}, please check your permissions!"
                )
            iresult, _ = result
            # do not use {} syntax below, it gives rise to a dict if the schema modules
            # are empty and then triggers a TypeError in missing_members = source_module - target_module
            source_module = set(  # noqa
                [mod for mod in iresult["schema_str"].split(",") if mod != ""]
            )
            target_module = ln_setup.settings.instance.modules
            if not source_module.issubset(target_module):
                missing_members = source_module - target_module
                logger.warning(
                    f"source modules has additional modules: {missing_members}\nconsider mounting these registry modules to transfer all metadata"
                )
            cache_filepath.write_text(f"{iresult['lnid']}\n{iresult['schema_str']}")  # type: ignore
            settings_file = instance_settings_file(name, owner)
            db = update_db_using_local(iresult, settings_file)
        else:
            isettings = load_instance_settings(settings_file)
            db = isettings.db
            cache_filepath.write_text(f"{isettings.uid}\n{','.join(isettings.modules)}")  # type: ignore
        add_db_connection(db, instance)
        return QuerySet(model=cls, using=instance)

    def __get_module_name__(cls) -> str:
        schema_module_name = cls.__module__.split(".")[0]
        module_name = schema_module_name.replace("lnschema_", "")
        if module_name == "lamindb":
            module_name = "core"
        return module_name

    @deprecated("__get_module_name__")
    def __get_schema_name__(cls) -> str:
        return cls.__get_module_name__()

    def __get_name_with_module__(cls) -> str:
        module_name = cls.__get_module_name__()
        if module_name == "core":
            module_prefix = ""
        else:
            module_prefix = f"{module_name}."
        return f"{module_prefix}{cls.__name__}"

    @deprecated("__get_name_with_module__")
    def __get_name_with_schema__(cls) -> str:
        return cls.__get_name_with_module__()


class BasicRecord(models.Model, metaclass=Registry):
    """Basic metadata record.

    It has the same methods as Record, but doesn't have the additional fields.

    It's mainly used for LinkORMs and similar.
    """

    class Meta:
        abstract = True

    def __init__(self, *args, **kwargs):
        skip_validation = kwargs.pop("_skip_validation", False)
        if not args:
            if (
                issubclass(self.__class__, Record)
                and not self.__class__.__name__ == "Storage"
            ):
                from lamindb import context as run_context

                if run_context.space is not None:
                    kwargs["space"] = run_context.space
            if skip_validation:
                super().__init__(**kwargs)
            else:
                from ..core._settings import settings
                from .can_curate import CanCurate
                from .collection import Collection
                from .schema import Schema
                from .transform import Transform

                validate_fields(self, kwargs)

                # do not search for names if an id is passed; this is important
                # e.g. when synching ids from the notebook store to lamindb
                has_consciously_provided_uid = False
                if "_has_consciously_provided_uid" in kwargs:
                    has_consciously_provided_uid = kwargs.pop(
                        "_has_consciously_provided_uid"
                    )
                if (
                    isinstance(self, (CanCurate, Collection, Transform))
                    and settings.creation.search_names
                    and not has_consciously_provided_uid
                ):
                    name_field = getattr(self, "_name_field", "name")
                    exact_match = suggest_records_with_similar_names(
                        self, name_field, kwargs
                    )
                    if exact_match is not None:
                        if "version" in kwargs:
                            if kwargs["version"] is not None:
                                version_comment = " and version"
                                existing_record = self.__class__.filter(
                                    **{
                                        name_field: kwargs[name_field],
                                        "version": kwargs["version"],
                                    }
                                ).one_or_none()
                            else:
                                # for a versioned record, an exact name match is not a criterion
                                # for retrieving a record in case `version` isn't passed -
                                # we'd always pull out many records with exactly the same name
                                existing_record = None
                        else:
                            version_comment = ""
                            existing_record = exact_match
                        if existing_record is not None:
                            logger.important(
                                f"returning existing {self.__class__.__name__} record with same"
                                f" {name_field}{version_comment}: '{kwargs[name_field]}'"
                            )
                            if isinstance(self, Schema):
                                if existing_record.hash != kwargs["hash"]:
                                    raise ValueError(
                                        f"Schema name is already in use by schema with uid '{existing_record.uid}', please choose a different name."
                                    )
                            init_self_from_db(self, existing_record)
                            update_attributes(self, kwargs)
                            return None
                super().__init__(**kwargs)
                if isinstance(self, ValidateFields):
                    # this will trigger validation against django validators
                    try:
                        if hasattr(self, "clean_fields"):
                            self.clean_fields()
                        else:
                            self._Model__clean_fields()
                    except DjangoValidationError as e:
                        message = _format_django_validation_error(self, e)
                        raise FieldValidationError(message) from e
        elif len(args) != len(self._meta.concrete_fields):
            raise FieldValidationError(
                f"Use keyword arguments instead of positional arguments, e.g.: {self.__class__.__name__}(name='...')."
            )
        else:
            super().__init__(*args)
        track_current_key_and_name_values(self)

    def save(self, *args, **kwargs) -> Record:
        """Save.

        Always saves to the default database.
        """
        using_key = None
        if "using" in kwargs:
            using_key = kwargs["using"]
        db = self._state.db
        pk_on_db = self.pk
        artifacts: list = []
        if self.__class__.__name__ == "Collection" and self.id is not None:
            # when creating a new collection without being able to access artifacts
            artifacts = self.ordered_artifacts.list()
        pre_existing_record = None
        # consider records that are being transferred from other databases
        transfer_logs: dict[str, list[str]] = {
            "mapped": [],
            "transferred": [],
            "run": None,
        }
        if db is not None and db != "default" and using_key is None:
            if isinstance(self, IsVersioned):
                if not self.is_latest:
                    raise NotImplementedError(
                        "You are attempting to transfer a record that's not the latest in its version history. This is currently not supported."
                    )
            pre_existing_record = transfer_to_default_db(
                self, using_key, transfer_logs=transfer_logs
            )
        self._revises: IsVersioned
        if pre_existing_record is not None:
            init_self_from_db(self, pre_existing_record)
        else:
            check_key_change(self)
            check_name_change(self)
            try:
                # save versioned record in presence of self._revises
                if isinstance(self, IsVersioned) and self._revises is not None:
                    assert self._revises.is_latest  # noqa: S101
                    revises = self._revises
                    revises.is_latest = False
                    with transaction.atomic():
                        revises._revises = None  # ensure we don't start a recursion
                        revises.save()
                        super().save(*args, **kwargs)  # type: ignore
                    self._revises = None
                # save unversioned record
                else:
                    super().save(*args, **kwargs)
            except IntegrityError as e:
                error_msg = str(e)
                # two possible error messages for hash duplication
                # "duplicate key value violates unique constraint"
                # "UNIQUE constraint failed"
                if (
                    "UNIQUE constraint failed" in error_msg
                    or "duplicate key value violates unique constraint" in error_msg
                ) and "hash" in error_msg:
                    pre_existing_record = self.__class__.get(hash=self.hash)
                    logger.warning(
                        f"returning {self.__class__.__name__.lower()} with same hash: {pre_existing_record}"
                    )
                    init_self_from_db(self, pre_existing_record)
                else:
                    raise
            # call the below in case a user makes more updates to the record
            track_current_key_and_name_values(self)
        # perform transfer of many-to-many fields
        # only supported for Artifact and Collection records
        if db is not None and db != "default" and using_key is None:
            if self.__class__.__name__ == "Collection":
                if len(artifacts) > 0:
                    logger.info("transfer artifacts")
                    for artifact in artifacts:
                        artifact.save()
                    self.artifacts.add(*artifacts)
            if hasattr(self, "labels"):
                from copy import copy

                from lamindb.models._feature_manager import FeatureManager

                # here we go back to original record on the source database
                self_on_db = copy(self)
                self_on_db._state.db = db
                self_on_db.pk = pk_on_db  # manually set the primary key
                self_on_db.features = FeatureManager(self_on_db)  # type: ignore
                self.features._add_from(self_on_db, transfer_logs=transfer_logs)
                self.labels.add_from(self_on_db, transfer_logs=transfer_logs)
            for k, v in transfer_logs.items():
                if k != "run":
                    logger.important(f"{k} records: {', '.join(v)}")

        if (
            self.__class__.__name__
            in {
                "Artifact",
                "Transform",
                "Run",
                "ULabel",
                "Feature",
                "Schema",
                "Collection",
                "Reference",
            }
            and self._branch_code >= 1
        ):
            import lamindb as ln

            if ln.context.project is not None:
                self.projects.add(ln.context.project)
        return self

    def delete(self) -> None:
        """Delete."""
        # note that the logic below does not fire if a record is moved to the trash
        # the idea is that moving a record to the trash should move its entire version family
        # to the trash, whereas permanently deleting should default to only deleting a single record
        # of a version family
        # we can consider making it easy to permanently delete entire version families as well,
        # but that's for another time
        if isinstance(self, IsVersioned) and self.is_latest:
            new_latest = (
                self.__class__.objects.using(self._state.db)
                .filter(is_latest=False, uid__startswith=self.stem_uid)
                .order_by("-created_at")
                .first()
            )
            if new_latest is not None:
                new_latest.is_latest = True
                with transaction.atomic():
                    new_latest.save()
                    super().delete()  # type: ignore
                logger.warning(f"new latest version is {new_latest}")
                return None
        super().delete()


class Space(BasicRecord):
    """Spaces.

    You can use spaces to restrict access to records within an instance.

    All data in this registry is synced from `lamin.ai` to enable re-using spaces across instances.
    There is no need to manually create records.
    """

    id: int = models.SmallAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    name: str = models.CharField(max_length=100, db_index=True)
    """Name of space."""
    uid: str = CharField(
        editable=False,
        unique=True,
        max_length=12,
        default="00000000",
        db_default="00000000",
        db_index=True,
    )
    """Universal id."""
    description: str | None = CharField(null=True)
    """Description of space."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of run."""

    @overload
    def __init__(
        self,
        name: str,
        description: str | None = None,
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
        super().__init__(*args, **kwargs)


@doc_args(RECORD_REGISTRY_EXAMPLE)
class Record(BasicRecord, metaclass=Registry):
    """Metadata record.

    Every `Record` is a data model that comes with a registry in form of a SQL
    table in your database.

    Sub-classing `Record` creates a new registry while instantiating a `Record`
    creates a new record.

    {}

    `Record`'s metaclass is :class:`~lamindb.models.Registry`.

    `Record` inherits from Django's `Model` class. Why does LaminDB call it `Record`
    and not `Model`? The term `Record` can't lead to confusion with statistical,
    machine learning or biological models.
    """

    _branch_code: int = models.SmallIntegerField(db_index=True, default=1, db_default=1)
    """Whether record is on a branch, in archive or in trash.

    This dictates whether a record appears in queries & searches.

    Coding is as follows:

    - 3: template (hidden in queries & searches)
    - 2: draft (hidden in queries & searches)
    - 1: default (visible in queries & searches)
    - 0: archive (hidden, meant to be kept)
    - -1: trash (hidden, scheduled for deletion)

    Any integer higher than >3 codes a branch that's involved in a pull request.
    """
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1)
    """The space in which the record lives."""
    _aux: dict[str, Any] | None = JSONField(default=None, db_default=None, null=True)
    """Auxiliary field for dictionary-like metadata."""

    class Meta:
        abstract = True


def _format_django_validation_error(record: Record, e: DjangoValidationError):
    """Pretty print Django validation errors."""
    errors = {}
    if hasattr(e, "error_dict"):
        error_dict = e.error_dict
    else:
        error_dict = {"__all__": e.error_list}

    for field_name, error_list in error_dict.items():
        for error in error_list:
            if hasattr(error, "message"):
                msg = error.message
            else:
                msg = str(error)

            if field_name == "__all__":
                errors[field_name] = f"{colors.yellow(msg)}"
            else:
                current_value = getattr(record, field_name, None)
                errors[field_name] = (
                    f"{field_name}: {colors.yellow(current_value)} is not valid\n    → {msg}"
                )

    if errors:
        message = "\n  "
        for _, error in errors.items():
            message += error + "\n  "

        return message


def _get_record_kwargs(record_class) -> list[tuple[str, str]]:
    """Gets the parameters of a Record from the overloaded signature.

    Example:
        >>> get_record_params(bt.Organism)
        >>> [('name', 'str'), ('taxon_id', 'str | None'), ('scientific_name', 'str | None')]
    """
    source = inspect.getsource(record_class)

    # Find first overload that's not *db_args
    pattern = r"@overload\s+def __init__\s*\(([\s\S]*?)\):\s*\.{3}"
    overloads = re.finditer(pattern, source)

    for single_overload in overloads:
        params_block = single_overload.group(1)
        # This is an additional safety measure if the overloaded signature that we're
        # looking for is not at the top but a "db_args" constructor
        if "*db_args" in params_block:
            continue

        params = []
        for line in params_block.split("\n"):
            line = line.strip()
            if not line or "self" in line:
                continue

            # Extract name and type annotation
            # The regex pattern finds parameter definitions like:
            # Simple: name: str
            # With default: age: int = 0
            # With complex types: items: List[str] = []
            param_pattern = (
                r"(\w+)"  # Parameter name
                r"\s*:\s*"  # Colon with optional whitespace
                r"((?:[^=,]|"  # Type hint: either non-equals/comma chars
                r"(?<=\[)[^[\]]*"  # or contents within square brackets
                r"(?=\]))+)"  # looking ahead for closing bracket
                r"(?:\s*=\s*"  # Optional default value part
                r"([^,]+))?"  # Default value: anything but comma
            )
            match = re.match(param_pattern, line)
            if not match:
                continue

            name, type_str = match.group(1), match.group(2).strip()

            # Keep type as string instead of evaluating
            params.append((name, type_str))

        return params

    return []


def _search(
    cls,
    string: str,
    *,
    field: StrField | list[StrField] | None = None,
    limit: int | None = 20,
    case_sensitive: bool = False,
    truncate_string: bool = False,
) -> QuerySet:
    if string is None:
        raise ValueError("Cannot search for None value! Please pass a valid string.")

    input_queryset = (
        cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    )
    registry = input_queryset.model
    name_field = getattr(registry, "_name_field", "name")
    if field is None:
        fields = [
            field.name
            for field in registry._meta.fields
            if field.get_internal_type() in {"CharField", "TextField"}
        ]
    else:
        if not isinstance(field, list):
            fields_input = [field]
        else:
            fields_input = field
        fields = []
        for field in fields_input:
            if not isinstance(field, str):
                try:
                    fields.append(field.field.name)
                except AttributeError as error:
                    raise TypeError(
                        "Please pass a Record string field, e.g., `CellType.name`!"
                    ) from error
            else:
                fields.append(field)

    if truncate_string:
        if (len_string := len(string)) > 5:
            n_80_pct = int(len_string * 0.8)
            string = string[:n_80_pct]

    string = string.strip()
    string_escape = re.escape(string)

    exact_lookup = Exact if case_sensitive else IExact
    regex_lookup = Regex if case_sensitive else IRegex
    contains_lookup = Contains if case_sensitive else IContains

    ranks = []
    contains_filters = []
    for field in fields:
        field_expr = Coalesce(
            Cast(field, output_field=TextField()),
            Value(""),
            output_field=TextField(),
        )
        # exact rank
        exact_expr = exact_lookup(field_expr, string)
        exact_rank = Cast(exact_expr, output_field=IntegerField()) * 200
        ranks.append(exact_rank)
        # exact synonym
        synonym_expr = regex_lookup(field_expr, rf"(?:^|.*\|){string_escape}(?:\|.*|$)")
        synonym_rank = Cast(synonym_expr, output_field=IntegerField()) * 200
        ranks.append(synonym_rank)
        # match as sub-phrase
        sub_expr = regex_lookup(
            field_expr, rf"(?:^|.*[ \|\.,;:]){string_escape}(?:[ \|\.,;:].*|$)"
        )
        sub_rank = Cast(sub_expr, output_field=IntegerField()) * 10
        ranks.append(sub_rank)
        # startswith and avoid matching string with " " on the right
        # mostly for truncated
        startswith_expr = regex_lookup(
            field_expr, rf"(?:^|.*\|){string_escape}[^ ]*(?:\|.*|$)"
        )
        startswith_rank = Cast(startswith_expr, output_field=IntegerField()) * 8
        ranks.append(startswith_rank)
        # match as sub-phrase from the left, mostly for truncated
        right_expr = regex_lookup(field_expr, rf"(?:^|.*[ \|]){string_escape}.*")
        right_rank = Cast(right_expr, output_field=IntegerField()) * 2
        ranks.append(right_rank)
        # match as sub-phrase from the right
        left_expr = regex_lookup(field_expr, rf".*{string_escape}(?:$|[ \|\.,;:].*)")
        left_rank = Cast(left_expr, output_field=IntegerField()) * 2
        ranks.append(left_rank)
        # simple contains filter
        contains_expr = contains_lookup(field_expr, string)
        contains_filter = Q(contains_expr)
        contains_filters.append(contains_filter)
        # also rank by contains
        contains_rank = Cast(contains_expr, output_field=IntegerField())
        ranks.append(contains_rank)
        # additional rule for truncated strings
        # weight matches from the beginning of the string higher
        # sometimes whole words get truncated and startswith_expr is not enough
        if truncate_string and field == name_field:
            startswith_lookup = StartsWith if case_sensitive else IStartsWith
            name_startswith_expr = startswith_lookup(field_expr, string)
            name_startswith_rank = (
                Cast(name_startswith_expr, output_field=IntegerField()) * 2
            )
            ranks.append(name_startswith_rank)

    ranked_queryset = (
        input_queryset.filter(reduce(lambda a, b: a | b, contains_filters))
        .alias(rank=sum(ranks))
        .order_by("-rank")
    )

    return ranked_queryset[:limit]


def _lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
    using_key: str | None = None,
) -> NamedTuple:
    """{}"""  # noqa: D415
    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    field = get_name_field(registry=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
    ).lookup(
        return_field=(
            get_name_field(registry=queryset.model, field=return_field)
            if return_field is not None
            else None
        )
    )


def get_name_field(
    registry: type[Record] | QuerySet | Manager,
    *,
    field: StrField | None = None,
) -> str:
    """Get the 1st char or text field from the registry."""
    if isinstance(registry, (QuerySet, Manager)):
        registry = registry.model
    model_field_names = [i.name for i in registry._meta.fields]

    # set to default name field
    if field is None:
        if hasattr(registry, "_name_field"):
            field = registry._meta.get_field(registry._name_field)
        elif "name" in model_field_names:
            field = registry._meta.get_field("name")
        else:
            # first char or text field that doesn't contain "id"
            for i in registry._meta.fields:
                if "id" in i.name:
                    continue
                if i.get_internal_type() in {"CharField", "TextField"}:
                    field = i
                    break

        # no default name field can be found
        if field is None:
            raise ValueError(
                "please pass a Record string field, e.g., `CellType.name`!"
            )
        else:
            field = field.name  # type:ignore
    if not isinstance(field, str):
        try:
            field = field.field.name
        except AttributeError:
            raise TypeError(
                "please pass a Record string field, e.g., `CellType.name`!"
            ) from None

    return field


def add_db_connection(db: str, using: str):
    db_config = dj_database_url.config(
        default=db, conn_max_age=600, conn_health_checks=True
    )
    db_config["TIME_ZONE"] = "UTC"
    db_config["OPTIONS"] = {}
    db_config["AUTOCOMMIT"] = True
    connections.settings[using] = db_config


REGISTRY_UNIQUE_FIELD = {
    "storage": "root",
    "feature": "name",
    "ulabel": "name",
    "space": "name",  # TODO: this should be updated with the currently used space instead during transfer
}


def update_fk_to_default_db(
    records: Record | list[Record] | QuerySet,
    fk: str,
    using_key: str | None,
    transfer_logs: dict,
):
    record = records[0] if isinstance(records, (list, QuerySet)) else records
    if hasattr(record, f"{fk}_id") and getattr(record, f"{fk}_id") is not None:
        fk_record = getattr(record, fk)
        field = REGISTRY_UNIQUE_FIELD.get(fk, "uid")
        fk_record_default = fk_record.__class__.filter(
            **{field: getattr(fk_record, field)}
        ).one_or_none()
        if fk_record_default is None:
            from copy import copy

            fk_record_default = copy(fk_record)
            transfer_to_default_db(
                fk_record_default, using_key, save=True, transfer_logs=transfer_logs
            )
        if isinstance(records, (list, QuerySet)):
            for r in records:
                setattr(r, f"{fk}", None)
                setattr(r, f"{fk}_id", fk_record_default.id)
        else:
            setattr(records, f"{fk}", None)
            setattr(records, f"{fk}_id", fk_record_default.id)


FKBULK = [
    "organism",
    "source",
    "report",  # Run
]


def transfer_fk_to_default_db_bulk(
    records: list | QuerySet, using_key: str | None, transfer_logs: dict
):
    for fk in FKBULK:
        update_fk_to_default_db(records, fk, using_key, transfer_logs=transfer_logs)


def get_transfer_run(record) -> Run:
    from lamindb import settings
    from lamindb.core._context import context
    from lamindb.models import Run, Transform
    from lamindb.models.artifact import WARNING_RUN_TRANSFORM

    slug = record._state.db
    owner, name = get_owner_name_from_identifier(slug)
    cache_filepath = ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    if not cache_filepath.exists():
        raise SystemExit("Need to call .using() before")
    instance_uid = cache_filepath.read_text().split("\n")[0]
    key = f"transfers/{instance_uid}"
    uid = instance_uid + "0000"
    transform = Transform.filter(uid=uid).one_or_none()
    if transform is None:
        search_names = settings.creation.search_names
        settings.creation.search_names = False
        transform = Transform(  # type: ignore
            uid=uid, description=f"Transfer from `{slug}`", key=key, type="function"
        ).save()
        settings.creation.search_names = search_names
    # use the global run context to get the initiated_by_run run id
    if context.run is not None:
        initiated_by_run = context.run
    else:
        if not settings.creation.artifact_silence_missing_run_warning:
            logger.warning(WARNING_RUN_TRANSFORM)
        initiated_by_run = None
    # it doesn't seem to make sense to create new runs for every transfer
    run = Run.filter(
        transform=transform, initiated_by_run=initiated_by_run
    ).one_or_none()
    if run is None:
        run = Run(transform=transform, initiated_by_run=initiated_by_run).save()  # type: ignore
        run.initiated_by_run = initiated_by_run  # so that it's available in memory
    return run


def transfer_to_default_db(
    record: Record,
    using_key: str | None,
    *,
    transfer_logs: dict,
    save: bool = False,
    transfer_fk: bool = True,
) -> Record | None:
    if record._state.db is None or record._state.db == "default":
        return None
    registry = record.__class__
    record_on_default = registry.objects.filter(uid=record.uid).one_or_none()
    record_str = f"{record.__class__.__name__}(uid='{record.uid}')"
    if transfer_logs["run"] is None:
        transfer_logs["run"] = get_transfer_run(record)
    if record_on_default is not None:
        transfer_logs["mapped"].append(record_str)
        return record_on_default
    else:
        transfer_logs["transferred"].append(record_str)

    if hasattr(record, "created_by_id"):
        record.created_by = None
        record.created_by_id = ln_setup.settings.user.id
    # run & transform
    run = transfer_logs["run"]
    if hasattr(record, "run_id"):
        record.run = None
        record.run_id = run.id
    # deal with denormalized transform FK on artifact and collection
    if hasattr(record, "transform_id"):
        record.transform = None
        record.transform_id = run.transform_id
    # transfer other foreign key fields
    fk_fields = [
        i.name
        for i in record._meta.fields
        if i.get_internal_type() == "ForeignKey"
        if i.name not in {"created_by", "run", "transform"}
    ]
    if not transfer_fk:
        # don't transfer fk fields that are already bulk transferred
        fk_fields = [fk for fk in fk_fields if fk not in FKBULK]
    for fk in fk_fields:
        update_fk_to_default_db(record, fk, using_key, transfer_logs=transfer_logs)
    record.id = None
    record._state.db = "default"
    if save:
        record.save()
    return None


def track_current_key_and_name_values(record: Record):
    from lamindb.models import Artifact

    if isinstance(record, Artifact):
        record._old_key = record.key
        record._old_suffix = record.suffix
    elif hasattr(record, "_name_field"):
        record._old_name = getattr(record, record._name_field)


def check_name_change(record: Record):
    """Warns if a record's name has changed."""
    from lamindb.models import Artifact, Collection, Feature, Schema, Transform

    if (
        not record.pk
        or not hasattr(record, "_old_name")
        or not hasattr(record, "_name_field")
    ):
        return

    # checked in check_key_change or not checked at all
    if isinstance(record, (Artifact, Collection, Transform)):
        return

    # renaming feature sets is not checked
    if isinstance(record, Schema):
        return

    old_name = record._old_name
    new_name = getattr(record, record._name_field)
    registry = record.__class__.__name__

    if old_name != new_name:
        # when a label is renamed, only raise a warning if it has a feature
        if hasattr(record, "artifacts"):
            linked_records = (
                record.artifacts.through.filter(
                    label_ref_is_name=True, **{f"{registry.lower()}_id": record.pk}
                )
                .exclude(feature_id=None)  # must have a feature
                .exclude(
                    feature_ref_is_name=None
                )  # must be linked via Curator and therefore part of a schema
                .distinct()
            )
            artifact_ids = linked_records.list("artifact__uid")
            n = len(artifact_ids)
            if n > 0:
                s = "s" if n > 1 else ""
                logger.error(
                    f"You are trying to {colors.red('rename label')} from '{old_name}' to '{new_name}'!\n"
                    f"   → The following {n} artifact{s} {colors.red('will no longer be validated')}: {artifact_ids}\n\n"
                    f"{colors.bold('To rename this label')}, make it external:\n"
                    f"   → run `artifact.labels.make_external(label)`\n\n"
                    f"After renaming, consider re-curating the above artifact{s}:\n"
                    f'   → in each dataset, manually modify label "{old_name}" to "{new_name}"\n'
                    f"   → run `ln.Curator`\n"
                )
                raise RecordNameChangeIntegrityError

        # when a feature is renamed
        elif isinstance(record, Feature):
            # only internal features are associated with schemas
            linked_artifacts = Artifact.filter(feature_sets__features=record).list(
                "uid"
            )
            n = len(linked_artifacts)
            if n > 0:
                s = "s" if n > 1 else ""
                logger.error(
                    f"You are trying to {colors.red('rename feature')} from '{old_name}' to '{new_name}'!\n"
                    f"   → The following {n} artifact{s} {colors.red('will no longer be validated')}: {linked_artifacts}\n\n"
                    f"{colors.bold('To rename this feature')}, make it external:\n"
                    "   → run `artifact.features.make_external(feature)`\n\n"
                    f"After renaming, consider re-curating the above artifact{s}:\n"
                    f"   → in each dataset, manually modify feature '{old_name}' to '{new_name}'\n"
                    f"   → run `ln.Curator`\n"
                )
                raise RecordNameChangeIntegrityError


def check_key_change(record: Union[Artifact, Transform]):
    """Errors if a record's key has falsely changed."""
    from .artifact import Artifact

    if not isinstance(record, Artifact) or not hasattr(record, "_old_key"):
        return
    if record._old_suffix != record.suffix:
        raise InvalidArgument(
            f"Changing the `.suffix` of an artifact is not allowed! You tried to change it from '{record._old_suffix}' to '{record.suffix}'."
        )

    old_key = record._old_key
    new_key = record.key

    if old_key != new_key:
        if not record._key_is_virtual:
            raise InvalidArgument(
                f"Changing a non-virtual key of an artifact is not allowed! You tried to change it from '{old_key}' to '{new_key}'."
            )
        if old_key is not None:
            old_key_suffix = extract_suffix_from_path(
                PurePosixPath(old_key), arg_name="key"
            )
            assert old_key_suffix == record.suffix, (  # noqa: S101
                old_key_suffix,
                record.suffix,
            )
        else:
            old_key_suffix = record.suffix
        new_key_suffix = extract_suffix_from_path(
            PurePosixPath(new_key), arg_name="key"
        )
        if old_key_suffix != new_key_suffix:
            raise InvalidArgument(
                f"The suffix '{new_key_suffix}' of the provided key is incorrect, it should be '{old_key_suffix}'."
            )


def format_field_value(value: datetime | str | Any) -> Any:
    from datetime import datetime

    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M:%S %Z")

    if isinstance(value, str):
        try:
            value = datetime.fromisoformat(value)
            value = value.strftime("%Y-%m-%d %H:%M:%S %Z")
        except ValueError:
            pass
        return f"'{value}'"
    else:
        return value


class RecordInfo:
    def __init__(self, registry: Registry):
        self.registry = registry

    def _get_type_for_field(self, field_name: str) -> str:
        field = self.registry._meta.get_field(field_name)
        related_model_name = (
            field.related_model.__name__
            if hasattr(field, "related_model") and field.related_model
            else None
        )
        return related_model_name if related_model_name else field.get_internal_type()

    def _get_base_class_fields(self) -> list[str]:
        return [
            field.name
            for base in self.registry.__bases__
            if hasattr(base, "_meta")
            for field in base._meta.get_fields()
        ]

    def _reorder_fields_by_class(self, fields_to_order: list[Field]) -> list[Field]:
        """Reorders the fields so that base class fields come last."""
        non_base_class_fields = [
            field
            for field in fields_to_order
            if field.name not in self._get_base_class_fields()
        ]
        found_base_class_fields = [
            field
            for field in fields_to_order
            if field.name in self._get_base_class_fields()
        ]
        return non_base_class_fields + found_base_class_fields

    def get_simple_fields(self, return_str: bool = False) -> Any:
        simple_fields = [
            field
            for field in self.registry._meta.get_fields()
            if not (
                isinstance(field, ManyToOneRel)
                or isinstance(field, ManyToManyRel)
                or isinstance(field, ManyToManyField)
                or isinstance(field, ForeignKey)
                or field.name.startswith("_")
                or field.name == "id"
            )
        ]
        simple_fields = self._reorder_fields_by_class(simple_fields)
        if not return_str:
            return simple_fields
        else:
            repr_str = f"  {colors.italic('Simple fields')}\n"
            if simple_fields:
                repr_str += "".join(
                    [
                        f"    .{field_name.name}: {self._get_type_for_field(field_name.name)}\n"
                        for field_name in simple_fields
                    ]
                )
            return repr_str

    def get_relational_fields(self, return_str: bool = False):
        # we ignore ManyToOneRel because it leads to so much clutter in the API
        # also note that our general guideline is to have related_name="+"
        # for ForeignKey fields
        relational_fields = (ManyToOneRel, ManyToManyRel, ManyToManyField, ForeignKey)

        class_specific_relational_fields = [
            field
            for field in self.registry._meta.fields + self.registry._meta.many_to_many
            if isinstance(field, relational_fields)
            and not field.name.startswith(("links_", "_"))
        ]

        non_class_specific_relational_fields = [
            field
            for field in self.registry._meta.get_fields()
            if isinstance(field, relational_fields)
            and not field.name.startswith(("links_", "_"))
        ]
        non_class_specific_relational_fields = self._reorder_fields_by_class(
            non_class_specific_relational_fields
        )

        # Ensure that class specific fields (e.g. Artifact) come before non-class specific fields (e.g. collection)
        filtered_non_class_specific = [
            field
            for field in non_class_specific_relational_fields
            if field not in class_specific_relational_fields
        ]
        ordered_relational_fields = (
            class_specific_relational_fields + filtered_non_class_specific
        )

        core_module_fields = []
        external_modules_fields = []
        for field in ordered_relational_fields:
            field_name = repr(field).split(": ")[1][:-1]
            if field_name.count(".") == 1 and "lamindb" not in field_name:
                external_modules_fields.append(field)
            else:
                core_module_fields.append(field)

        def _get_related_field_type(field) -> str:
            field_type = (
                field.related_model.__get_name_with_module__()
                .replace(
                    "Artifact", ""
                )  # some fields have an unnecessary 'Artifact' in their name
                .replace(
                    "Collection", ""
                )  # some fields have an unnecessary 'Collection' in their name
            )
            return (
                self._get_type_for_field(field.name)
                if not field_type.strip()
                else field_type
            )

        core_module_fields_formatted = [
            f"    .{field.name}: {_get_related_field_type(field)}\n"
            for field in core_module_fields
        ]
        external_modules_fields_formatted = [
            f"    .{field.name}: {_get_related_field_type(field)}\n"
            for field in external_modules_fields
        ]

        if not return_str:
            external_modules_fields_by_modules = defaultdict(list)
            for field_str, field in zip(
                external_modules_fields_formatted, external_modules_fields
            ):
                field_type = field_str.split(":")[1].split()[0]
                module_name = field_type.split(".")[0]
                external_modules_fields_by_modules[module_name].append(field)
            return core_module_fields, external_modules_fields_by_modules
        else:
            repr_str = ""

            # Non-external relational fields
            if core_module_fields:
                repr_str += f"  {colors.italic('Relational fields')}\n"
                repr_str += "".join(core_module_fields_formatted)

            # External relational fields
            external_modules = set()
            for field in external_modules_fields_formatted:
                field_type = field.split(":")[1].split()[0]
                external_modules.add(field_type.split(".")[0])

            if external_modules:
                # We want Bionty to show up before other modules
                external_modules = (
                    ["bionty"] + sorted(external_modules - {"bionty"})  # type: ignore
                    if "bionty" in external_modules
                    else sorted(external_modules)
                )
                for ext_module in external_modules:
                    ext_module_fields = [
                        field
                        for field in external_modules_fields_formatted
                        if ext_module in field
                    ]

                    if ext_module_fields:
                        repr_str += (
                            f"  {colors.italic(f'{ext_module.capitalize()} fields')}\n"
                        )
                        repr_str += "".join(ext_module_fields)

            return repr_str


def registry_repr(cls):
    """Shows fields."""
    repr_str = f"{colors.green(cls.__name__)}\n"
    info = RecordInfo(cls)
    repr_str += info.get_simple_fields(return_str=True)
    repr_str += info.get_relational_fields(return_str=True)
    repr_str = repr_str.rstrip("\n")
    return repr_str


def record_repr(
    self: Record, include_foreign_keys: bool = True, exclude_field_names=None
) -> str:
    if exclude_field_names is None:
        exclude_field_names = ["id", "updated_at", "source_code"]
    field_names = [
        field.name
        for field in self._meta.fields
        if (not isinstance(field, ForeignKey) and field.name not in exclude_field_names)
    ]
    if include_foreign_keys:
        field_names += [
            f"{field.name}_id"
            for field in self._meta.fields
            if isinstance(field, ForeignKey)
        ]
    if "created_at" in field_names:
        field_names.remove("created_at")
        field_names.append("created_at")
    if field_names[0] != "uid" and "uid" in field_names:
        field_names.remove("uid")
        field_names.insert(0, "uid")
    fields_str = {}
    for k in field_names:
        if not k.startswith("_") and hasattr(self, k):
            value = getattr(self, k)
            # Force strip the time component of the version
            if k == "version" and value:
                fields_str[k] = f"'{str(value).split()[0]}'"
            else:
                fields_str[k] = format_field_value(value)
    fields_joined_str = ", ".join(
        [f"{k}={fields_str[k]}" for k in fields_str if fields_str[k] is not None]
    )
    return f"{self.__class__.__name__}({fields_joined_str})"


# below is code to further format the repr of a record
#
# def format_repr(
#     record: Record, exclude_field_names: str | list[str] | None = None
# ) -> str:
#     if isinstance(exclude_field_names, str):
#         exclude_field_names = [exclude_field_names]
#     exclude_field_names_init = ["id", "created_at", "updated_at"]
#     if exclude_field_names is not None:
#         exclude_field_names_init += exclude_field_names
#     return record.__repr__(
#         include_foreign_keys=False, exclude_field_names=exclude_field_names_init
#     )


Record.__repr__ = record_repr  # type: ignore
Record.__str__ = record_repr  # type: ignore


class Migration(BasicRecord):
    app = CharField(max_length=255)
    name = CharField(max_length=255)
    applied: datetime = DateTimeField()

    class Meta:
        db_table = "django_migrations"
        managed = False
