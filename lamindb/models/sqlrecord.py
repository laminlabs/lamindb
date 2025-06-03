from __future__ import annotations

import builtins
import inspect
import re
import sys
from collections import defaultdict
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
from django.db import IntegrityError, ProgrammingError, connections, models, transaction
from django.db.models import CASCADE, PROTECT, Field, Manager, QuerySet
from django.db.models.base import ModelBase
from django.db.models.fields.related import (
    ManyToManyField,
    ManyToManyRel,
    ManyToOneRel,
)
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
    update_db_using_local,
)
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance_hub
from lamindb_setup.core._settings_store import instance_settings_file
from lamindb_setup.core.django import DBToken, db_token_manager
from lamindb_setup.core.upath import extract_suffix_from_path

from lamindb.base import deprecated

from ..base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
)
from ..base.ids import base62_12
from ..base.types import FieldAttr, StrField
from ..errors import (
    FieldValidationError,
    InvalidArgument,
    NoWriteAccess,
    SQLRecordNameChangeIntegrityError,
    ValidationError,
)
from ._is_versioned import IsVersioned
from .query_manager import QueryManager, _lookup, _search

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd

    from .artifact import Artifact
    from .run import Run, User
    from .transform import Transform


T = TypeVar("T", bound="SQLRecord")
IPYTHON = getattr(builtins, "__IPYTHON__", False)


# -------------------------------------------------------------------------------------
# A note on required fields at the SQLRecord level
#
# As Django does most of its validation on the Form-level, it doesn't offer functionality
# for validating the integrity of an SQLRecord object upon instantation (similar to pydantic)
#
# For required fields, we define them as commonly done on the SQL level together
# with a validator in SQLRecord (validate_required_fields)
#
# This goes against the Django convention, but goes with the SQLModel convention
# (Optional fields can be null on the SQL level, non-optional fields cannot)
#
# Due to Django's convention where CharFieldAttr has pre-configured (null=False, default=""), marking
# a required field necessitates passing `default=None`. Without the validator it would trigger
# an error at the SQL-level, with it, it triggers it at instantiation

# -------------------------------------------------------------------------------------
# A note on class and instance methods of core SQLRecord
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


class IsLink:
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


def init_self_from_db(self: SQLRecord, existing_record: SQLRecord):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def update_attributes(record: SQLRecord, attributes: dict[str, str]):
    for key, value in attributes.items():
        if getattr(record, key) != value and value is not None:
            if key not in {"uid", "dtype", "otype", "hash"}:
                logger.warning(f"updated {key} from {getattr(record, key)} to {value}")
                setattr(record, key, value)
            else:
                hash_message = (
                    "recomputing on .save()"
                    if key == "hash"
                    else f"keeping {getattr(record, key)}"
                )
                logger.warning(
                    f"ignoring tentative value {value} for {key}, {hash_message}"
                )


def validate_literal_fields(record: SQLRecord, kwargs) -> None:
    """Validate all Literal type fields in a record.

    Args:
        record: record being validated

    Raises:
        ValidationError: If any field value is not in its Literal's allowed values
    """
    if isinstance(record, IsLink):
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


def validate_fields(record: SQLRecord, kwargs):
    from lamindb.models import (
        Artifact,
        Collection,
        Feature,
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
    }:
        uid_max_length = record.__class__._meta.get_field(
            "uid"
        ).max_length  # triggers FieldDoesNotExist
        if len(kwargs["uid"]) != uid_max_length:  # triggers KeyError
            if not (
                record.__class__ is Schema and len(kwargs["uid"]) == 16
            ):  # no error for schema
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
    record: SQLRecord, name_field: str, kwargs
) -> SQLRecord | None:
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

        from lamindb import SQLRecord, fields

        # sub-classing `SQLRecord` creates a new registry
        class Experiment(SQLRecord):
            name: str = fields.CharField()

        # instantiating `Experiment` creates a record `experiment`
        experiment = Experiment(name="my experiment")

        # you can save the record to the database
        experiment.save()

        # `Experiment` refers to the registry, which you can query
        df = Experiment.filter(name__startswith="my ").df()
"""


# this is the metaclass for SQLRecord
@doc_args(RECORD_REGISTRY_EXAMPLE)
class Registry(ModelBase):
    """Metaclass for :class:`~lamindb.models.SQLRecord`.

    Each `Registry` *object* is a `SQLRecord` *class* and corresponds to a table in the metadata SQL database.

    You work with `Registry` objects whenever you use *class methods* of `SQLRecord`.

    You call any subclass of `SQLRecord` a "registry" and their objects "records". A `SQLRecord` object corresponds to a row in the SQL table.

    If you want to create a new registry, you sub-class `SQLRecord`.

    {}

    Note: `Registry` inherits from Django's `ModelBase`.
    """

    _available_fields: set[str] = None

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

    @doc_args(_lookup.__doc__)
    def lookup(
        cls,
        field: StrField | None = None,
        return_field: StrField | None = None,
    ) -> NamedTuple:
        """{}"""  # noqa: D415
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

        Examples:

            ::

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

    @doc_args(_search.__doc__)
    def search(
        cls,
        string: str,
        *,
        field: StrField | None = None,
        limit: int | None = 20,
        case_sensitive: bool = False,
    ) -> QuerySet:
        """{}"""  # noqa: D415
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

        # connection already established
        if instance in connections:
            return QuerySet(model=cls, using=instance)
        # we're in the default instance
        if instance is None or instance == "default":
            return QuerySet(model=cls, using=None)
        owner, name = get_owner_name_from_identifier(instance)
        if [owner, name] == setup_settings.instance.slug.split("/"):
            return QuerySet(model=cls, using=None)

        # move on to different instances
        cache_using_filepath = (
            setup_settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
        )
        settings_file = instance_settings_file(name, owner)
        if not settings_file.exists():
            result = connect_instance_hub(owner=owner, name=name)
            if isinstance(result, str):
                raise RuntimeError(
                    f"Failed to load instance {instance}, please check your permissions!"
                )
            iresult, _ = result
            # do not use {} syntax below, it gives rise to a dict if the schema modules
            # are empty and then triggers a TypeError in missing_members = source_modules - target_modules
            source_modules = set(  # noqa
                [mod for mod in iresult["schema_str"].split(",") if mod != ""]
            )
            # this just retrives the full connection string from iresult
            db = update_db_using_local(iresult, settings_file)
            cache_using_filepath.write_text(
                f"{iresult['lnid']}\n{iresult['schema_str']}"
            )
            # need to set the token if it is a fine_grained_access and the user is jwt (not public)
            is_fine_grained_access = (
                iresult["fine_grained_access"] and iresult["db_permissions"] == "jwt"
            )
            # access_db can take both: the dict from connect_instance_hub and isettings
            into_db_token = iresult
        else:
            isettings = load_instance_settings(settings_file)
            source_modules = isettings.modules
            db = isettings.db
            cache_using_filepath.write_text(
                f"{isettings.uid}\n{','.join(source_modules)}"
            )
            # need to set the token if it is a fine_grained_access and the user is jwt (not public)
            is_fine_grained_access = (
                isettings._fine_grained_access and isettings._db_permissions == "jwt"
            )
            # access_db can take both: the dict from connect_instance_hub and isettings
            into_db_token = isettings

        target_modules = setup_settings.instance.modules
        if missing_members := source_modules - target_modules:
            logger.info(
                f"in transfer, source lamindb instance has additional modules: {', '.join(missing_members)}"
            )

        add_db_connection(db, instance)
        if is_fine_grained_access:
            db_token = DBToken(into_db_token)
            db_token_manager.set(db_token, instance)
        return QuerySet(model=cls, using=instance)

    def __get_module_name__(cls) -> str:
        schema_module_name = cls.__module__.split(".")[0]
        module_name = schema_module_name.replace("lnschema_", "")
        if module_name == "lamindb":
            module_name = "core"
        return module_name

    def __get_name_with_module__(cls) -> str:
        module_name = cls.__get_module_name__()
        if module_name == "core":
            module_prefix = ""
        else:
            module_prefix = f"{module_name}."
        return f"{module_prefix}{cls.__name__}"

    def __get_available_fields__(cls) -> set[str]:
        if cls._available_fields is None:
            cls._available_fields = {
                f.name
                for f in cls._meta.get_fields()
                if not f.name.startswith("_")
                and not f.name.startswith("links_")
                and not f.name.endswith("_id")
            }
            if cls.__name__ == "Artifact":
                cls._available_fields.add("visibility")  # backward compat
                cls._available_fields.add("_branch_code")  # backward compat
                cls._available_fields.add("transform")
        return cls._available_fields


class BaseSQLRecord(models.Model, metaclass=Registry):
    """Basic metadata record.

    It has the same methods as SQLRecord, but doesn't have the additional fields.

    It's mainly used for IsLinks and similar.
    """

    objects = QueryManager()

    class Meta:
        abstract = True
        base_manager_name = "objects"

    def __init__(self, *args, **kwargs):
        skip_validation = kwargs.pop("_skip_validation", False)
        if not args:
            if self.__class__.__name__ in {
                "Artifact",
                "Collection",
                "Transform",
                "Run",
            }:
                from lamindb import context as run_context

                if run_context.space is not None:
                    kwargs["space"] = run_context.space
            if issubclass(
                self.__class__, SQLRecord
            ) and self.__class__.__name__ not in {"Storage", "Source"}:
                from lamindb import context as run_context

                if run_context.branch is not None:
                    kwargs["branch"] = run_context.branch
            if skip_validation:
                super().__init__(**kwargs)
            else:
                from ..core._settings import settings
                from .can_curate import CanCurate
                from .collection import Collection
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

    def save(self, *args, **kwargs) -> SQLRecord:
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
            except (IntegrityError, ProgrammingError) as e:
                error_msg = str(e)
                # two possible error messages for hash duplication
                # "duplicate key value violates unique constraint"
                # "UNIQUE constraint failed"
                if (
                    self.__class__.__name__ in {"Transform", "Artifact"}
                    and isinstance(e, IntegrityError)
                    and "hash" in error_msg
                    and (
                        "UNIQUE constraint failed" in error_msg
                        or "duplicate key value violates unique constraint" in error_msg
                    )
                ):
                    pre_existing_record = self.__class__.get(hash=self.hash)
                    logger.warning(
                        f"returning {self.__class__.__name__.lower()} with same hash: {pre_existing_record}"
                    )
                    init_self_from_db(self, pre_existing_record)
                elif (
                    isinstance(e, ProgrammingError)
                    and hasattr(self, "space")
                    and "new row violates row-level security policy" in error_msg
                ):
                    raise NoWriteAccess(
                        f"You’re not allowed to write to the space '{self.space.name}'.\n"
                        "Please contact an administrator of the space if you need write access."
                    ) from None
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
                if k != "run" and len(v) > 0:
                    logger.important(f"{k} records: {', '.join(v)}")

        if self.__class__.__name__ in {
            "Artifact",
            "Transform",
            "Run",
            "ULabel",
            "Feature",
            "Schema",
            "Collection",
            "Reference",
        } and not (
            self.__class__.__name__ == "Artifact" and self.kind == "__lamindb_run__"
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


class Space(BaseSQLRecord):
    """Workspaces with managed access for specific users or teams.

    Guide: :doc:`docs:access`.

    All data in this registry is synchronized from LaminHub so that spaces can be shared and reused across multiple LaminDB instances.
    """

    id: int = models.SmallAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    name: str = models.CharField(max_length=100, db_index=True)
    """Name of space."""
    uid: str = CharField(
        editable=False,
        unique=True,
        max_length=12,
        default="A",
        db_default="A",
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
    """Creator of space."""

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


class Branch(BaseSQLRecord):
    """Branches for change management with archive and trash states.

    Every `SQLRecord` has a `branch` field, which dictates where a record appears in queries & searches.
    """

    # below isn't fully implemented but a roadmap
    # - 3: template (hidden in queries & searches)
    # - 2: locked (same as default, but locked for edits except for space admins)
    # - 1: default (visible in queries & searches)
    # - 0: archive (hidden, meant to be kept, locked for edits for everyone)
    # - -1: trash (hidden, scheduled for deletion)

    # An integer higher than >3 codes a branch that can be used for collaborators to create drafts
    # that can be merged onto the main branch in an experience akin to a Pull Request. The mapping
    # onto a semantic branch name is handled through LaminHub.

    id: int = models.AutoField(primary_key=True)
    """An integer id that's synchronized for a family of coupled database instances.

    Among all LaminDB instances, this id is arbitrary and non-unique.
    """
    name: str = models.CharField(max_length=100, db_index=True)
    """Name of branch."""
    uid: str = CharField(
        editable=False,
        unique=True,
        max_length=12,
        default=base62_12,
        db_index=True,
    )
    """Universal id.

    This id is useful if one wants to apply the same patch to many database instances.
    """
    description: str | None = CharField(null=True)
    """Description of branch."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of branch."""

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
class SQLRecord(BaseSQLRecord, metaclass=Registry):
    """Metadata record.

    Every `SQLRecord` is a data model that comes with a registry in form of a SQL
    table in your database.

    Sub-classing `SQLRecord` creates a new registry while instantiating a `SQLRecord`
    creates a new record.

    {}

    `SQLRecord`'s metaclass is :class:`~lamindb.models.Registry`.

    `SQLRecord` inherits from Django's `Model` class. Why does LaminDB call it `SQLRecord`
    and not `Model`? The term `SQLRecord` can't lead to confusion with statistical,
    machine learning or biological models.
    """

    branch: Branch = ForeignKey(
        Branch,
        PROTECT,
        default=1,
        db_default=1,
        db_column="_branch_code",
        related_name="+",
    )
    """Whether record is on a branch or in another "special state"."""
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1, related_name="+")
    """The space in which the record lives."""
    _aux: dict[str, Any] | None = JSONField(default=None, db_default=None, null=True)
    """Auxiliary field for dictionary-like metadata."""

    class Meta:
        abstract = True

    @property
    @deprecated("branch_id")
    def _branch_code(self) -> int:
        """Deprecated alias for `branch`."""
        return self.branch_id

    @_branch_code.setter
    def _branch_code(self, value: int):
        self.branch_id = value


def _format_django_validation_error(record: SQLRecord, e: DjangoValidationError):
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
    """Gets the parameters of a SQLRecord from the overloaded signature.

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


def get_name_field(
    registry: type[SQLRecord] | QuerySet | Manager,
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
                "please pass a SQLRecord string field, e.g., `CellType.name`!"
            )
        else:
            field = field.name  # type:ignore
    if not isinstance(field, str):
        try:
            field = field.field.name
        except AttributeError:
            raise TypeError(
                "please pass a SQLRecord string field, e.g., `CellType.name`!"
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


REGISTRY_UNIQUE_FIELD = {"storage": "root", "feature": "name", "ulabel": "name"}


def update_fk_to_default_db(
    records: SQLRecord | list[SQLRecord] | QuerySet,
    fk: str,
    using_key: str | None,
    transfer_logs: dict,
):
    # here in case it is an iterable, we are checking only a single record
    # and set the same fks for all other records because we do this only
    # for certain fks where they have to the same for the whole bulk
    # see transfer_fk_to_default_db_bulk
    # todo: but this has to be changed i think, it is not safe as it is now - Sergei
    record = records[0] if isinstance(records, (list, QuerySet)) else records
    if getattr(record, f"{fk}_id", None) is not None:
        # set the space of the transferred record to the current space
        if fk == "space":
            # for space we set the record's space to the current space
            from lamindb import context

            # the default space has id=1
            fk_record_default = Space.get(1) if context.space is None else context.space
        # process non-space fks
        else:
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
        # re-set the fks to the newly saved ones in the default db
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
    cache_using_filepath = (
        ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    )
    if not cache_using_filepath.exists():
        raise SystemExit("Need to call .using() before")
    instance_uid = cache_using_filepath.read_text().split("\n")[0]
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
    record: SQLRecord,
    using_key: str | None,
    *,
    transfer_logs: dict,
    save: bool = False,
    transfer_fk: bool = True,
) -> SQLRecord | None:
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
        if i.name not in {"created_by", "run", "transform", "branch"}
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


def track_current_key_and_name_values(record: SQLRecord):
    from lamindb.models import Artifact

    # below, we're using __dict__ to avoid triggering the refresh from the database
    # which can lead to a recursion
    if isinstance(record, Artifact):
        record._old_key = record.__dict__.get("key")
        record._old_suffix = record.__dict__.get("suffix")
    elif hasattr(record, "_name_field"):
        record._old_name = record.__dict__.get(record._name_field)


def check_name_change(record: SQLRecord):
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
                raise SQLRecordNameChangeIntegrityError

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
                raise SQLRecordNameChangeIntegrityError


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


class SQLRecordInfo:
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
    info = SQLRecordInfo(cls)
    repr_str += info.get_simple_fields(return_str=True)
    repr_str += info.get_relational_fields(return_str=True)
    repr_str = repr_str.rstrip("\n")
    return repr_str


def record_repr(
    self: SQLRecord, include_foreign_keys: bool = True, exclude_field_names=None
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
#     record: SQLRecord, exclude_field_names: str | list[str] | None = None
# ) -> str:
#     if isinstance(exclude_field_names, str):
#         exclude_field_names = [exclude_field_names]
#     exclude_field_names_init = ["id", "created_at", "updated_at"]
#     if exclude_field_names is not None:
#         exclude_field_names_init += exclude_field_names
#     return record.__repr__(
#         include_foreign_keys=False, exclude_field_names=exclude_field_names_init
#     )


SQLRecord.__repr__ = record_repr  # type: ignore
SQLRecord.__str__ = record_repr  # type: ignore


class Migration(BaseSQLRecord):
    app = CharField(max_length=255)
    name = CharField(max_length=255)
    applied: datetime = DateTimeField()

    class Meta:
        db_table = "django_migrations"
        managed = False


LinkORM = IsLink  # backward compat
Record = SQLRecord  # backward compat
BasicRecord = BaseSQLRecord  # backward compat
RecordInfo = SQLRecordInfo  # backward compat
