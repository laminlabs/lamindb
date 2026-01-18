from __future__ import annotations

import builtins
import gzip
import inspect
import os
import re
import shutil
import sys
from collections import defaultdict
from itertools import chain
from pathlib import Path, PurePosixPath
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
from django.db.models import ForeignKey as django_ForeignKey
from django.db.models.base import ModelBase
from django.db.models.fields.related import (
    ManyToManyField,
    ManyToManyRel,
    ManyToOneRel,
)
from django.db.models.functions import Lower
from lamin_utils import colors, logger
from lamindb_setup import settings as setup_settings
from lamindb_setup._connect_instance import (
    INSTANCE_NOT_FOUND_MESSAGE,
    InstanceNotFoundError,
    get_owner_name_from_identifier,
    load_instance_settings,
    update_db_using_local,
)
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance_hub
from lamindb_setup.core._settings_store import instance_settings_file
from lamindb_setup.core.django import DBToken, db_token_manager
from lamindb_setup.core.upath import extract_suffix_from_path
from upath import UPath

from lamindb.base.utils import class_and_instance_method, deprecated

from ..base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    TextField,
)
from ..base.types import FieldAttr, StrField
from ..base.uids import base62_12
from ..errors import (
    FieldValidationError,
    InvalidArgument,
    NoWriteAccess,
    ValidationError,
)
from ._is_versioned import IsVersioned
from .query_manager import QueryManager, _lookup, _search

if TYPE_CHECKING:
    from datetime import datetime

    import pandas as pd

    from .artifact import Artifact
    from .block import BranchBlock, SpaceBlock
    from .query_set import SQLRecordList
    from .run import Run, User
    from .transform import Transform


T = TypeVar("T", bound="SQLRecord")
IPYTHON = getattr(builtins, "__IPYTHON__", False)
UNIQUE_FIELD_NAMES = {
    "root",
    "ontology_id",
    "uid",
    "scientific_name",
    "ensembl_gene_id",
    "uniprotkb_id",
}


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


class IsLink:
    pass


class HasType(models.Model):
    """Mixin for registries that have a hierarchical `type` assigned.

    Such registries have a `.type` foreign key pointing to themselves.

    A `type` hence allows hierarchically grouping records under types.

    For instance, using the example of `ln.Record`::

        experiment_type = ln.Record(name="Experiment", is_type=True).save()
        experiment1 = ln.Record(name="Experiment 1", type=experiment_type).save()
        experiment2 = ln.Record(name="Experiment 2", type=experiment_type).save()
    """

    class Meta:
        abstract = True

    is_type: bool = BooleanField(default=False, db_default=False, db_index=True)
    """Indicates if record is a `type`.

    For example, if a record "Compound" is a `type`, the actual compounds "darerinib", "tramerinib", would be instances of that `type`.
    """

    def query_types(self) -> SQLRecordList:
        """Query types of a record recursively.

        While `.type` retrieves the `type`, this method
        retrieves all super types of that `type`::

            # Create type hierarchy
            type1 = model_class(name="Type1", is_type=True).save()
            type2 = model_class(name="Type2", is_type=True, type=type1).save()
            type3 = model_class(name="Type3", is_type=True, type=type2).save()

            # Create a record with type3
            record = model_class(name=f"{model_name}3", type=type3).save()

            # Query super types
            super_types = record.query_types()
            assert super_types[0] == type3
            assert super_types[1] == type2
            assert super_types[2] == type1
        """
        from .has_parents import _query_ancestors_of_fk

        return _query_ancestors_of_fk(self, "type")  # type: ignore


def deferred_attribute__repr__(self):
    return f"FieldAttr({self.field.model.__name__}.{self.field.name})"


def unique_constraint_error_in_error_message(error_msg: str) -> bool:
    """Check if the error message indicates a unique constraint violation."""
    return (
        "UNIQUE constraint failed" in error_msg  # SQLite
        or "duplicate key value violates unique constraint" in error_msg  # Postgre
    )


def parse_violated_field_from_error_message(error_msg: str) -> list[str] | None:
    # Even if the model has multiple fields with unique=True,
    # Django will only raise an IntegrityError for one field at a time
    # - whichever constraint is violated first during the database insert/update operation.
    if unique_constraint_error_in_error_message(error_msg):
        if "UNIQUE constraint failed" in error_msg:  # sqlite
            constraint_field = (
                error_msg.removeprefix("UNIQUE constraint failed: ")
                .split(", ")[0]
                .split(".")[-1]
            )
            return [constraint_field]
        else:  # postgres
            # Extract constraint name from double quotes
            constraint_name = error_msg.split('"')[1]

            # Check if it's a multi-column constraint (contains multiple field names)
            # Format: tablename_field1_field2_..._hash_uniq
            if "_uniq" in constraint_name:
                # Remove '_uniq' suffix first
                constraint_name = constraint_name.removesuffix("_uniq")

                # Remove hash (8 hex characters at the end)
                parts = constraint_name.split("_")
                if len(parts[-1]) == 8 and all(
                    c in "0123456789abcdef" for c in parts[-1]
                ):
                    constraint_name = "_".join(parts[:-1])

                # Remove table name prefix (e.g., "bionty_ethnicity_")
                # Table name is typically the first 2 parts for app_model format
                parts = constraint_name.split("_")
                if len(parts) > 2:
                    # Assume first 2 parts are table name (e.g., "bionty_ethnicity")
                    field_string = "_".join(parts[2:])
                else:
                    field_string = constraint_name

                # Now parse the fields from DETAIL line
                # DETAIL: Key (name, ontology_id)=(South Asian, HANCESTRO:0006) already exists.
                if "Key (" in error_msg:
                    fields_part = error_msg.split("Key (")[1].split(")=")[0]
                    fields = [f.strip() for f in fields_part.split(",")]
                    return fields

                # Fallback if DETAIL line not available
                return [field_string]
            else:
                # Single field constraint (ends with _key)
                constraint_field = constraint_name.removesuffix("_key").split("_")[-1]
                return [constraint_field]

    return None


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
    from .run import current_run

    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"
    # if run was not set on the existing record, set it to the current_run
    if hasattr(self, "run_id") and self.run_id is None and current_run() is not None:
        logger.warning(f"run was not set on {self}, setting to current run")
        self.run = current_run()


def update_attributes(record: SQLRecord, attributes: dict[str, str]):
    for key, value in attributes.items():
        if getattr(record, key) != value and value is not None:
            if key not in {"uid", "_dtype_str", "otype", "hash"}:
                logger.warning(f"updated {key} from {getattr(record, key)} to {value}")
                setattr(record, key, value)
            else:
                hash_message = (
                    "recomputing on .save()"
                    if key == "hash"
                    else f"keeping {getattr(record, key)}"
                )
                logger.debug(
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
    from lamindb.base.types import ArtifactKind, Dtype, TransformKind

    types = {
        "TransformKind": TransformKind,
        "ArtifactKind": ArtifactKind,
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
    if (
        "type" in kwargs
        and isinstance(kwargs["type"], HasType)
        and not kwargs["type"].is_type
    ):
        object_name = record.__class__.__name__.lower()
        raise ValueError(
            f"You can only assign a {object_name} with `is_type=True` as `type` to another {object_name}, but this doesn't have it: {kwargs['type']}"
        )
    # validate literals
    validate_literal_fields(record, kwargs)


def suggest_records_with_similar_names(
    record: SQLRecord, name_field: str, kwargs
) -> SQLRecord | None:
    """Returns a record if found exact match, otherwise None.

    Logs similar matches if found.
    """
    if kwargs.get(name_field) is None or not isinstance(kwargs.get(name_field), str):
        return None
    # need to perform an additional request to find the exact match
    # previously, this was inferred from the truncated/fuzzy search below
    # but this isn't reliable: https://laminlabs.slack.com/archives/C04FPE8V01W/p1737812808563409
    # the below needs to be .first() because there might be multiple records with the same
    # name field in case the record is versioned (e.g. for Transform key)
    if isinstance(record, HasType):
        if kwargs.get("type", None) is None:
            subset = record.__class__.filter(type__isnull=True)
        else:
            subset = record.__class__.filter(type=kwargs["type"])
    else:
        subset = record.__class__
    exact_match = subset.filter(**{name_field: kwargs[name_field]}).first()
    if exact_match is not None:
        return exact_match
    queryset = _search(
        subset,
        kwargs[name_field],
        field=name_field,
        truncate_string=True,
        limit=3,
    )
    if not queryset.exists():  # empty queryset
        return None
    s, it, nots, record_text = (
        ("", "it", "s", "a record")
        if len(queryset) == 1
        else ("s", "one of them", "", "records")
    )
    similar_names = ", ".join(f"'{getattr(record, name_field)}'" for record in queryset)
    msg = f"you are trying to create a record with name='{kwargs[name_field]}' but {record_text} with similar {name_field}{s} exist{nots}: {similar_names}. Did you mean to load {it}?"
    logger.warning(f"{msg}")

    return None


def delete_record(record: BaseSQLRecord, is_soft: bool = True):
    def delete():
        if is_soft:
            record.branch_id = -1
            record.save()
            return None
        else:
            return super(BaseSQLRecord, record).delete()

    # deal with versioned records
    # if _ovewrite_version = True, there is only a single version and
    # no need to set the new latest version because all versions are deleted
    # when deleting the latest version
    if (
        isinstance(record, IsVersioned)
        and record.is_latest
        and not getattr(record, "_overwrite_versions", False)
    ):
        new_latest = (
            record.__class__.objects.using(record._state.db)
            .filter(is_latest=False, uid__startswith=record.stem_uid)
            .exclude(branch_id=-1)  # exclude candidates in the trash
            .order_by("-created_at")
            .first()
        )
        if new_latest is not None:
            new_latest.is_latest = True
            if is_soft:
                record.is_latest = False
            with transaction.atomic():
                new_latest.save()
                result = delete()
            logger.important_hint(f"new latest version is: {new_latest}")
            return result
    # deal with all other cases of the nested if condition now
    return delete()


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
        df = Experiment.filter(name__startswith="my ").to_dataframe()
"""


def _synchronize_clone(storage_root: str) -> str | None:
    """Synchronizes a clone to the local SQLite path.

    Args:
        storage_root: The storage root path of the (target) instance
    """
    cloud_db_path = UPath(storage_root) / ".lamindb" / "lamin.db"
    local_sqlite_path = ln_setup.settings.cache_dir / cloud_db_path.path.lstrip("/")

    if local_sqlite_path.exists():
        return f"sqlite:///{local_sqlite_path}"

    local_sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    cloud_db_path_gz = UPath(str(cloud_db_path) + ".gz", anon=True)
    local_sqlite_path_gz = Path(str(local_sqlite_path) + ".gz")

    try:
        cloud_db_path_gz.synchronize_to(
            local_sqlite_path_gz, error_no_origin=True, print_progress=True
        )
        with gzip.open(local_sqlite_path_gz, "rb") as f_in:
            with open(local_sqlite_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return f"sqlite:///{local_sqlite_path}"
    except (FileNotFoundError, PermissionError):
        logger.debug("Clone not found. Falling back to normal access...")
        return None


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

    def describe(cls, return_str: bool = False) -> str | None:
        """Describe the fields of the registry."""
        from ._describe import strip_ansi_from_string as _strip_ansi

        repr_str = f"{colors.green(cls.__name__)}\n"
        info = SQLRecordInfo(cls)
        repr_str += info.get_simple_fields(return_str=True)
        repr_str += info.get_relational_fields(return_str=True)
        repr_str = repr_str.rstrip("\n")
        if return_str:
            return _strip_ansi(repr_str)
        else:
            print(repr_str)
            return None

    @doc_args(_lookup.__doc__)
    def lookup(
        cls,
        field: StrField | None = None,
        return_field: StrField | None = None,
        keep: Literal["first", "last", False] = "first",
    ) -> NamedTuple:
        """{}"""  # noqa: D415
        return _lookup(cls=cls, field=field, return_field=return_field, keep=keep)

    def filter(cls, *queries, **expressions) -> QuerySet:
        """Query records.

        Args:
            queries: One or multiple `Q` objects.
            expressions: Fields and values passed as Django query expressions.

        See Also:
            - Guide: :doc:`docs:registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples:
            >>> ln.Project(name="my label").save()
            >>> ln.Project.filter(name__startswith="my").to_dataframe()
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
            :exc:`lamindb.errors.ObjectDoesNotExist`: In case no matching record is found.

        See Also:
            - Guide: :doc:`registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples:

            ::

                record = ln.Record.get("FvtpPJLJ")
                record = ln.Record.get(name="my-label")
        """
        from .query_set import QuerySet

        return QuerySet(model=cls).get(idlike, **expressions)

    def to_dataframe(
        cls,
        *,
        include: str | list[str] | None = None,
        features: str | list[str] | None = None,
        limit: int | None = 100,
        order_by: str | None = "-id",
    ) -> pd.DataFrame:
        """Evaluate and convert to `pd.DataFrame`.

        By default, maps simple fields and foreign keys onto `DataFrame` columns.

        Guide: :doc:`docs:registries`

        Args:
            include: Related data to include as columns. Takes strings of
                form `"records__name"`, `"cell_types__name"`, etc. or a list
                of such strings. For `Artifact`, `Record`, and `Run`, can also pass `"features"`
                to include features with data types pointing to entities in the core schema.
                If `"privates"`, includes private fields (fields starting with `_`).
            features: Configure the features to include. Can be a feature name or a list of such names.
                If `"queryset"`, infers the features used within the current queryset.
                Only available for `Artifact`, `Record`, and `Run`.
            limit: Maximum number of rows to display. If `None`, includes all results.
            order_by: Field name to order the records by. Prefix with '-' for descending order.
                Defaults to '-id' to get the most recent records. This argument is ignored
                if the queryset is already ordered or if the specified field does not exist.

        Examples:

            Include the name of the creator::

                ln.Record.to_dataframe(include="created_by__name"])

            Include features::

                ln.Artifact.to_dataframe(include="features")

            Include selected features::

                ln.Artifact.to_dataframe(features=["cell_type_by_expert", "cell_type_by_model"])
        """
        return cls.filter().to_dataframe(
            include=include, features=features, order_by=order_by, limit=limit
        )

    @deprecated(new_name="to_dataframe")
    def df(
        cls,
        *,
        include: str | list[str] | None = None,
        features: str | list[str] | None = None,
        limit: int | None = 100,
        order_by: str | None = "-id",
    ) -> pd.DataFrame:
        return cls.to_dataframe(
            include=include, features=features, limit=limit, order_by=order_by
        )

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

    @deprecated(new_name="connect")
    def using(
        cls,
        instance: str | None,
    ) -> QuerySet:
        return cls.connect(
            instance=instance,
        )

    def connect(
        cls,
        instance: str | None,
    ) -> QuerySet:
        """Query a non-default LaminDB instance.

        Args:
            instance: An instance identifier of form "account_handle/instance_name".

        Examples:

            ::

                ln.Record.connect("account_handle/instance_name").search("label7", field="name")
        """
        from .query_set import QuerySet

        # we're in the default instance
        if instance is None or instance == "default":
            return QuerySet(model=cls, using=None)
        # connection already established
        if instance in connections:
            return QuerySet(model=cls, using=instance)

        owner, name = get_owner_name_from_identifier(instance)
        current_instance_owner_name: list[str] = setup_settings.instance.slug.split("/")

        # move on to different instances
        cache_using_filepath = (
            setup_settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
        )
        settings_file = instance_settings_file(name, owner)
        if not settings_file.exists():
            result = connect_instance_hub(owner=owner, name=name)
            if isinstance(result, str):
                message = INSTANCE_NOT_FOUND_MESSAGE.format(
                    owner=owner, name=name, hub_result=result
                )
                raise InstanceNotFoundError(message)
            iresult, storage = result
            # this can happen if querying via an old instance name
            if [iresult.get("owner"), iresult["name"]] == current_instance_owner_name:
                return QuerySet(model=cls, using=None)
            # do not use {} syntax below, it gives rise to a dict if the schema modules
            # are empty and then triggers a TypeError in missing_members = source_modules - target_modules
            source_modules = set(  # noqa
                [mod for mod in iresult["schema_str"].split(",") if mod != ""]
            )

            # Try to connect to a clone if targeting a public instance but fall back to normal access if access failed
            db = None
            if (
                "_public" in iresult["db_user_name"]
                and "postgresql" in iresult["db_scheme"]
            ):
                db = _synchronize_clone(storage["root"])
            if db is None:
                if [
                    iresult.get("owner"),
                    iresult["name"],
                ] == current_instance_owner_name:
                    return QuerySet(model=cls, using=None)
                db = update_db_using_local(iresult, settings_file)
                is_fine_grained_access = (
                    iresult["fine_grained_access"]
                    and iresult["db_permissions"] == "jwt"
                )
            else:
                is_fine_grained_access = False

            cache_using_filepath.write_text(
                f"{iresult['lnid']}\n{iresult['schema_str']}", encoding="utf-8"
            )

            # access_db can take both: the dict from connect_instance_hub and isettings
            into_db_token = iresult
        else:
            isettings = load_instance_settings(settings_file)
            source_modules = isettings.modules
            db = None
            if "public" in isettings.db and isettings.dialect == "postgresql":
                db = _synchronize_clone(isettings.storage.root_as_str)

            # Try to connect to a clone if targeting a public instance but fall back to normal access if access failed
            if db is None:
                if [isettings.owner, isettings.name] == current_instance_owner_name:
                    return QuerySet(model=cls, using=None)
                db = isettings.db
                is_fine_grained_access = (
                    isettings._fine_grained_access
                    and isettings._db_permissions == "jwt"
                )
            else:
                is_fine_grained_access = False

            cache_using_filepath.write_text(
                f"{isettings.uid}\n{','.join(source_modules)}", encoding="utf-8"
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
            available_fields = set()
            for field in cls._meta.get_fields():
                if not (field_name := field.name).startswith(("_", "links_")):
                    available_fields.add(field_name)
                    if isinstance(field, django_ForeignKey):
                        available_fields.add(field_name + "_id")
            if cls.__name__ == "Artifact":
                available_fields.add("transform")
                available_fields.add("feature_sets")  # backward compat with lamindb v1
            cls._available_fields = available_fields
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

    # fields to track for changes
    # if not None, will be tracked in self._original_values as {field_name: value}
    # use _id fields for foreign keys
    _TRACK_FIELDS: tuple[str, ...] | None = None

    def __init__(self, *args, **kwargs):
        skip_validation = kwargs.pop("_skip_validation", False)
        if not args:
            if not os.getenv("LAMINDB_MULTI_INSTANCE") == "true":
                if (
                    issubclass(self.__class__, SQLRecord)
                    and self.__class__.__name__ != "Storage"
                    # do not save bionty entities in restricted spaces by default
                    and self.__class__.__module__ != "bionty.models"
                ):
                    from lamindb import context as run_context

                    if run_context.space is not None:
                        current_space = run_context.space
                    elif setup_settings.space is not None:
                        current_space = setup_settings.space

                    if current_space is not None:
                        if "space_id" in kwargs:
                            # space_id takes precedence over space
                            # https://claude.ai/share/f045e5dc-0143-4bc5-b8a4-38309229f75e
                            if kwargs["space_id"] == 1:  # ignore default space
                                kwargs.pop("space_id")
                                kwargs["space"] = current_space
                        elif "space" in kwargs:
                            if kwargs["space"] is None:
                                kwargs["space"] = current_space
                        else:
                            kwargs["space"] = current_space
                if issubclass(
                    self.__class__, SQLRecord
                ) and self.__class__.__name__ not in {"Storage", "Source"}:
                    from lamindb import context as run_context

                    if run_context.branch is not None:
                        current_branch = run_context.branch
                    elif setup_settings.branch is not None:
                        current_branch = setup_settings.branch

                    if current_branch is not None:
                        # branch_id takes precedence over branch
                        # https://claude.ai/share/f045e5dc-0143-4bc5-b8a4-38309229f75e
                        if "branch_id" in kwargs:
                            if kwargs["branch_id"] == 1:  # ignore default branch
                                kwargs.pop("branch_id")
                                kwargs["branch"] = current_branch
                        elif "branch" in kwargs:
                            if kwargs["branch"] is None:
                                kwargs["branch"] = current_branch
                        else:
                            kwargs["branch"] = current_branch
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
                        if "version_tag" in kwargs:
                            if kwargs.get("version_tag") is not None:
                                version_comment = " and version"
                                existing_record = self.__class__.filter(
                                    **{
                                        name_field: kwargs[name_field],
                                        "version_tag": kwargs.get("version_tag"),
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
                                f"returning {self.__class__.__name__.lower()} with same"
                                f" {name_field}{version_comment}: '{kwargs[name_field]}'"
                            )
                            init_self_from_db(self, existing_record)
                            update_attributes(self, kwargs)
                            # track original values after replacing with the existing record
                            self._populate_tracked_fields()
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
        # track original values of fields that are tracked for changes
        self._populate_tracked_fields()
        # TODO: refactor to use _TRACK_FIELDS
        track_current_key_and_name_values(self)

    # used in __init__
    # populates the _original_values dictionary with the original values of the tracked fields
    def _populate_tracked_fields(self):
        if (track_fields := self._TRACK_FIELDS) is not None:
            self._original_values = {f: self.__dict__[f] for f in track_fields}
        else:
            self._original_values = None

    def _field_changed(self, field_name: str) -> bool:
        """Check if the field has changed since the record was saved."""
        # use _id fields for foreign keys in field_name
        if self._state.adding:
            return False
        # check if the field is tracked for changes
        track_fields = self._TRACK_FIELDS
        assert track_fields is not None, (
            "_TRACK_FIELDS must be set for the record to track changes"
        )
        assert field_name in track_fields, (
            f"Field {field_name} is not tracked for changes"
        )
        # check if the field has changed since the record was created
        return self._original_values[field_name] != self.__dict__[field_name]

    def save(self: T, *args, **kwargs) -> T:
        """Save.

        Always saves to the default database.
        """
        using_key = None
        if "using" in kwargs:
            using_key = kwargs["using"]
        transfer_config = kwargs.pop("transfer", None)
        db = self._state.db
        pk_on_db = self.pk
        artifacts: list = []
        if self.__class__.__name__ == "Collection" and self.id is not None:
            # when creating a new collection without being able to access artifacts
            artifacts = self.ordered_artifacts.to_list()
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
            # TODO: refactor to use _TRACK_FIELDS
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
                # error for hash/uid duplication
                if (
                    self.__class__.__name__ in {"Transform", "Artifact", "Collection"}
                    and isinstance(e, IntegrityError)
                    and "hash" in error_msg
                    and unique_constraint_error_in_error_message(error_msg)
                ):
                    # we also need to include the key here because hash can be the same across keys
                    query_fields = {"hash": self.hash, "key": self.key}
                    if self.__class__.__name__ == "Artifact":
                        # in case of artifact, also storage is needed
                        query_fields["storage"] = self.storage
                    # the get here is Django's get and not aware of the trash or other branches
                    # but generally we bypass branch_id in queries for hash also in LaminDB's get()
                    pre_existing_record = self.__class__.get(**query_fields)
                    from_trash = (
                        "from trash" if pre_existing_record.branch_id == -1 else ""
                    )
                    pre_existing_record.branch_id = 1  # move to default branch
                    logger.warning(
                        f"returning {self.__class__.__name__.lower()} {from_trash} with same hash & key: {pre_existing_record}"
                    )
                    init_self_from_db(self, pre_existing_record)
                elif (
                    isinstance(e, IntegrityError)
                    # for Storage, even if uid was in the error message, we can retrieve based on
                    # the root because it's going to be the same root
                    and any(field in error_msg for field in UNIQUE_FIELD_NAMES)
                    and (
                        "_type_name_at_" not in error_msg
                    )  # constraints for unique type names in Record, ULabel, etc.
                    and (
                        "UNIQUE constraint failed" in error_msg
                        or "duplicate key value violates unique constraint" in error_msg
                    )
                    and hasattr(self, "branch_id")
                ):
                    unique_fields = parse_violated_field_from_error_message(error_msg)
                    # here we query against the all branches with .objects
                    pre_existing_record = self.__class__.objects.get(
                        **{field: getattr(self, field) for field in unique_fields}
                    )
                    # if the existing record is in the default branch, we just return it
                    if pre_existing_record.branch_id == 1:
                        logger.warning(
                            f"returning {self.__class__.__name__} record with same {unique_fields}: '{ {field: getattr(self, field) for field in unique_fields} }'"
                        )
                    # if the existing record is in a different branch we update its fields
                    else:
                        # modifies the fields of the existing record with new values of self
                        field_names = [i.name for i in self.__class__._meta.fields]
                        update_attributes(
                            pre_existing_record,
                            {f: getattr(self, f) for f in field_names},
                        )
                        pre_existing_record.save()
                    init_self_from_db(self, pre_existing_record)
                elif (
                    isinstance(e, ProgrammingError)
                    and "new row violates row-level security policy" in error_msg
                    and (
                        (is_locked := getattr(self, "is_locked", False))
                        or hasattr(self, "space")
                    )
                ):
                    if is_locked:
                        no_write_msg = "It is not allowed to modify or create locked ('is_locked=True') records."
                    else:
                        no_write_msg = (
                            f"You're not allowed to write to the space '{self.space.name}'.\n"
                            "Please contact administrators of the space if you need write access."
                        )
                    raise NoWriteAccess(no_write_msg) from None
                elif (
                    isinstance(e, ProgrammingError)
                    and "permission denied for table" in error_msg
                    and (isettings := setup_settings.instance)._db_permissions
                    == "public"
                ):
                    slug = isettings.slug
                    raise NoWriteAccess(
                        f"You are trying to write to '{slug}' with public (read-only) permissions.\n"
                        "Please contact administrators to make you a collaborator if you need write access.\n"
                        f"If you are already a collaborator, please do 'lamin connect {slug}' in console, "
                        "restart the python session and try again."
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
            if hasattr(self, "labels") and transfer_config == "annotations":
                from copy import copy

                # here we go back to original record on the source database
                self_on_db = copy(self)
                self_on_db._state.db = db
                self_on_db.pk = pk_on_db  # manually set the primary key
                self.features._add_from(self_on_db, transfer_logs=transfer_logs)
                self.labels.add_from(self_on_db, transfer_logs=transfer_logs)
            for k, v in transfer_logs.items():
                if k != "run" and len(v) > 0:
                    logger.important(f"{k}: {', '.join(v)}")

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

    @class_and_instance_method
    def describe(cls_or_self, return_str: bool = False) -> None | str:
        """Describe record including relations.

        Args:
            return_str: Return a string instead of printing.
        """
        from ._describe import describe_postgres_sqlite

        if isinstance(cls_or_self, type):
            return type(cls_or_self).describe(cls_or_self, return_str=return_str)  # type: ignore
        else:
            return describe_postgres_sqlite(cls_or_self, return_str=return_str)

    def __repr__(
        self: SQLRecord,
        include_foreign_keys: bool = True,
        exclude_field_names: list[str] | None = None,
    ) -> str:
        if exclude_field_names is None:
            exclude_field_names = ["id", "updated_at", "source_code"]
        field_names = [
            field.name
            for field in self._meta.fields
            if (
                not isinstance(field, ForeignKey)
                and field.name not in exclude_field_names
            )
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
        if "is_locked" in field_names:
            field_names.remove("is_locked")
            field_names.append("is_locked")
        if field_names[0] != "uid" and "uid" in field_names:
            field_names.remove("uid")
            field_names.insert(0, "uid")
        fields_str = {}
        for k in field_names:
            if k == "n" and getattr(self, k) < 0:
                # only needed for Schema
                continue
            if (
                not k.startswith("_")
                or (k == "_dtype_str" and self.__class__.__name__ == "Feature")
            ) and hasattr(self, k):
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

    def __str__(self) -> str:
        return self.__repr__()

    def delete(self, permanent: bool | None = None):
        """Delete.

        Args:
            permanent: For consistency, `False` raises an error, as soft delete is impossible.

        Returns:
            When `permanent=True`, returns Django's delete return value: a tuple of
            (deleted_count, {registry_name: count}). Otherwise returns None.
        """
        if permanent is False:
            raise ValueError(
                f"Soft delete is not possible for {self.__class__.__name__}, "
                "use 'permanent=True' or 'permanent=None' for permanent deletion."
            )

        return delete_record(self, is_soft=False)


class Space(BaseSQLRecord):
    """Spaces with managed access for specific users or teams.

    If not setting a space, a :class:`~lamindb.models.SQLRecord` object is accessible to all collaborators of the LaminDB instance because its :attr:`~lamindb.models.SQLRecord.space` field defaults to the built-in `all` space.
    You can create a restricted space through LaminHub either on the instance settings page or the *Spaces* tab of your account page.

    Examples:

        After creating a restricted space through LaminHub, create an artifact in the space::

            space = ln.Space.get(name="Our space")  # get a space
            ln.Artifact("./test.txt", key="test.txt", space=space).save()  # save artifact in space

        You can also move an existing object into a space::

            space = ln.Space.get(name="Our space")  # select a space
            record = ln.Record.get(name="existing label")
            record.space = space
            record.save()  # saved in space "Our space"

        For more examples and background, see :doc:`docs:permissions`, in particular, section :ref:`docs:use-a-restricted-space`.

    Notes:

        All data in this registry is synchronized from LaminHub so that spaces can be shared and reused across multiple LaminDB instances.
    """

    class Meta:
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(Lower("name"), name="unique_space_name_lower")
        ]

    id: int = models.SmallAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    name: str = models.CharField(max_length=100, db_index=True)
    """Name of space."""
    uid: str = CharField(
        editable=False,
        unique=True,
        max_length=12,
        default=base62_12,
        db_index=True,
    )
    """Universal id."""
    description: str | None = TextField(null=True)
    """Description of space."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of space."""
    ablocks: SpaceBlock
    """Blocks that annotate this space."""

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
        if not args and "uid" not in kwargs:
            warn = False
            msg = ""
            isettings = setup_settings.instance
            if (dialect := isettings.dialect) != "postgresql":
                warn = True
                msg = f"on {dialect} databases"
            elif not isettings.is_on_hub:
                warn = True
                msg = "on local instances"
            if warn:
                logger.warning(
                    f"creating spaces manually {msg} is possible for demo purposes, "
                    "but does *not* affect access permissions"
                )
        super().__init__(*args, **kwargs)


class Branch(BaseSQLRecord):
    """Branches for change management with archive and trash states.

    There are 3 pre-defined branches: `main`, `trash`, and `archive`.

    You can create branches similar to `git` via `lamin create --branch my_branch`.

    To add objects to that new branch rather than the `main` branch, run `lamin switch --branch my_branch`.

    To merge a set of artifacts on the `"my_branch"` branch to the main branch, run::

        ln.Artifact.filter(branch__name="my_branch").update(branch_id=1)

    If you delete an object via `sqlrecord.delete()`, it gets moved to the `trash` branch and scheduled for deletion.
    """

    class Meta:
        app_label = "lamindb"
        constraints = [
            models.UniqueConstraint(Lower("name"), name="unique_branch_name_lower")
        ]

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
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1, related_name="+")
    """The space associated with the branch."""
    description: str | None = TextField(null=True)
    """Description of branch."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of branch."""
    ablocks: BranchBlock
    """Blocks that annotate this branch."""

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
    """An object that maps to a row in a SQL table in the database.

    Every `SQLRecord` is a data model that comes with a registry in form of a SQL
    table in your database.

    Sub-classing `SQLRecord` creates a new registry while instantiating a `SQLRecord`
    creates a new object.

    {}

    `SQLRecord`'s metaclass is :class:`~lamindb.models.Registry`.

    `SQLRecord` inherits from Django's `Model` class. Why does LaminDB call it `SQLRecord`
    and not `Model`? The term `SQLRecord` can't lead to confusion with statistical,
    machine learning or biological models.
    """

    # we need the db_default when not interacting via django directly on a required field
    branch: Branch = ForeignKey(
        Branch,
        PROTECT,
        default=1,
        db_default=1,
        db_column="branch_id",
        related_name="+",
    )
    """The branch."""
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1, related_name="+")
    """The space."""
    is_locked: bool = BooleanField(default=False, db_default=False)
    """Whether the object is locked for edits."""
    _aux: dict[str, Any] | None = JSONField(default=None, db_default=None, null=True)
    """Auxiliary field for dictionary-like metadata."""

    class Meta:
        abstract = True

    def restore(self) -> None:
        """Restore from trash onto the main branch.

        Does **not** restore descendant objects if the object is `HasType` with `is_type = True`.
        """
        self.branch_id = 1
        self.save()

    def delete(self, permanent: bool | None = None, **kwargs):
        """Delete object.

        If object is `HasType` with `is_type = True`, deletes all descendant objects, too.

        Args:
            permanent: Whether to permanently delete the object (skips trash).
                If `None`, performs soft delete if the object is not already in the trash.

        Returns:
            When `permanent=True`, returns Django's delete return value: a tuple of
            (deleted_count, {registry_name: count}). Otherwise returns None.

        Examples:

            For any `SQLRecord` object `sqlrecord`, call::

                sqlrecord.delete()
        """
        if self._state.adding:
            logger.warning("record is not yet saved, delete has no effect")
            return None
        name_with_module = self.__class__.__get_name_with_module__()

        if name_with_module == "Artifact":
            # this first check means an invalid delete fails fast rather than cascading through
            # database and storage permission errors
            isettings = setup_settings.instance
            if self.storage.instance_uid != isettings.uid and (
                kwargs["storage"] or kwargs["storage"] is None
            ):
                from ..errors import IntegrityError
                from .storage import Storage

                raise IntegrityError(
                    "Cannot simply delete artifacts outside of this instance's managed storage locations."
                    "\n(1) If you only want to delete the metadata record in this instance, pass `storage=False`"
                    f"\n(2) If you want to delete the artifact in storage, please connect to the writing lamindb instance (uid={self.storage.instance_uid})."
                    f"\nThese are all managed storage locations of this instance:\n{Storage.filter(instance_uid=isettings.uid).to_dataframe()}"
                )

        # change branch_id to trash
        trash_branch_id = -1
        if self.branch_id > trash_branch_id and permanent is not True:
            if isinstance(self, HasType) and self.is_type:
                for child in getattr(
                    self, f"query_{self.__class__.__name__.lower()}s"
                )():
                    child.delete()
            delete_record(self, is_soft=True)
            logger.important(f"moved record to trash: {self}")
            return None

        # permanent delete
        if permanent is None:
            object_type_name = self.__class__.__name__
            log_identifier = self.uid if hasattr(self, "uid") else self.pk
            response = input(
                f"{object_type_name} {log_identifier} is already in trash! Are you sure you want to delete it from your"
                " database? You can't undo this action. (y/n) "
            )
            confirm_delete = response == "y"
        else:
            confirm_delete = permanent

        if confirm_delete:
            if name_with_module == "Run":
                from .run import delete_run_artifacts

                delete_run_artifacts(self)
            elif name_with_module == "Transform":
                from .transform import delete_transform_relations

                delete_transform_relations(self)
            elif name_with_module == "Artifact":
                from .artifact import delete_permanently

                delete_permanently(
                    self, storage=kwargs["storage"], using_key=kwargs["using_key"]
                )
                return None
            if name_with_module != "Artifact":
                return super().delete()
        return None


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
                f"Do not know which field to use as name file for registry {registry}, please pass field"
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


REGISTRY_UNIQUE_FIELD = {"storage": "root", "ulabel": "name"}


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
        raise SystemExit("Need to call .connect() before")
    instance_uid = cache_using_filepath.read_text().split("\n")[0]
    key = f"__lamindb_transfer__/{instance_uid}"
    uid = instance_uid + "0000"
    transform = Transform.filter(uid=uid).one_or_none()
    if transform is None:
        search_names = settings.creation.search_names
        settings.creation.search_names = False
        transform = Transform(  # type: ignore
            uid=uid, description=f"Transfer from `{slug}`", key=key, kind="function"
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
    run = Run.filter(transform=transform, initiated_by_run=initiated_by_run).first()
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
    logger.debug(f"transferring {registry.__name__} record {record.uid} to default db")
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
        record._old_key = record.__dict__.get("key")  # type: ignore
        record._old_suffix = record.__dict__.get("suffix")  # type: ignore
    elif hasattr(record, "_name_field"):
        record._old_name = record.__dict__.get(record._name_field)


def check_name_change(record: SQLRecord):
    """Warns if a record's name has changed."""
    from lamindb.models import (
        Artifact,
        Collection,
        Feature,
        Schema,
        Storage,
        Transform,
    )

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
        if hasattr(record, "artifacts") and not isinstance(record, Storage):
            linked_records = (
                # find all artifacts that are linked to this label via a feature with dtype
                # matching on the name aka "[registry]"
                record.artifacts.through.filter(
                    feature___dtype_str__contains=f"[{registry}]",
                    **{f"{registry.lower()}_id": record.pk},
                )
            )
            artifact_uids = list(set(linked_records.to_list("artifact__uid")))
            n = len(artifact_uids)
            if n > 0:
                s = "s" if n > 1 else ""
                es = "es" if n == 1 else ""
                logger.error(
                    f"by {colors.red('renaming label')} from '{old_name}' to '{new_name}' "
                    f"{n} artifact{s} no longer match{es} the label name in storage: {artifact_uids}\n\n"
                    f"   → consider re-curating\n"
                )
        elif isinstance(record, Feature):
            # only internal features of schemas with `itype=Feature` are prone to getting out of sync
            artifact_uids = Artifact.filter(
                schemas__features=record, schemas__itype="Feature"
            ).to_list("uid")
            n = len(artifact_uids)
            if n > 0:
                s = "s" if n > 1 else ""
                es = "es" if n == 1 else ""
                logger.warning(
                    f"by {colors.red('renaming feature')} from '{old_name}' to '{new_name}' "
                    f"{n} artifact{s} no longer match{es} the feature name in storage: {artifact_uids}\n"
                    "  → consider re-curating"
                )


def check_key_change(record: Union[Artifact, Transform]):
    """Errors if a record's key has falsely changed."""
    from .artifact import Artifact

    if not isinstance(record, Artifact) or not hasattr(record, "_old_key"):
        return
    if hasattr(record, "_skip_key_change_check") and record._skip_key_change_check:
        return
    if record._old_suffix != record.suffix:  # type: ignore
        raise InvalidArgument(
            f"Changing the `.suffix` of an artifact is not allowed! You tried to change it from '{record._old_suffix}' to '{record.suffix}'."  # type: ignore
        )

    old_key = record._old_key  # type: ignore
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


def format_field_value(value: datetime | str | Any, none: str = "None") -> str:
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
    if value is None:
        return none
    return str(value)


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

        # For Record class, move linked_in fields to the end
        if self.registry.__name__ == "Record":
            regular_fields = [
                f
                for f in ordered_relational_fields
                if not f.name.startswith(("linked_", "values_"))
            ]
            linked_fields = [
                f for f in ordered_relational_fields if f.name.startswith("linked_")
            ]
            values_fields = [
                f for f in ordered_relational_fields if f.name.startswith("values_")
            ]
            ordered_relational_fields = regular_fields + linked_fields + values_fields

        core_module_fields = []
        external_modules_fields = []
        for field in ordered_relational_fields:
            field_name = repr(field).split(": ")[1][:-1]
            if field_name.count(".") == 1 and "lamindb" not in field_name:
                external_modules_fields.append(field)
            else:
                core_module_fields.append(field)

        def _get_related_field_type(field) -> str:
            model_name = field.related_model.__get_name_with_module__()
            # Extract the class name (after the last dot if there's a module prefix)
            class_name = model_name.split(".")[-1]
            # Skip replacement for compound names like ArtifactBlock, FeatureBlock, etc.
            if class_name.endswith("Block"):
                # Return just the class name for Block types
                field_type = class_name
            else:
                field_type = (
                    model_name.replace(
                        "Artifact", ""
                    ).replace(  # some fields have an unnecessary 'Artifact' in their name
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


class Migration(BaseSQLRecord):
    app = CharField(max_length=255)
    name = CharField(max_length=255)
    applied: datetime = DateTimeField()

    class Meta:
        db_table = "django_migrations"
        app_label = "lamindb"
        managed = False


LinkORM = IsLink  # backward compat
Record = SQLRecord  # backward compat
BasicRecord = BaseSQLRecord  # backward compat
RecordInfo = SQLRecordInfo  # backward compat
