from __future__ import annotations

import builtins
import inspect
import re
from functools import reduce
from pathlib import PurePosixPath
from typing import TYPE_CHECKING, NamedTuple

import dj_database_url
import lamindb_setup as ln_setup
from django.core.exceptions import ValidationError as DjangoValidationError
from django.db import connections, transaction
from django.db.models import (
    IntegerField,
    Manager,
    Q,
    QuerySet,
    TextField,
    Value,
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
from django.db.utils import IntegrityError
from lamin_utils import colors, logger
from lamin_utils._lookup import Lookup
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
    update_db_using_local,
)
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance_hub
from lamindb_setup.core._settings_store import instance_settings_file
from lamindb_setup.core.upath import extract_suffix_from_path

from lamindb.errors import FieldValidationError
from lamindb.models import (
    Artifact,
    BasicRecord,
    CanCurate,
    Collection,
    Feature,
    IsVersioned,
    Param,
    Record,
    Run,
    Schema,
    Transform,
    ULabel,
    ValidateFields,
)

from ._utils import attach_func_to_class_method
from .core._settings import settings
from .errors import (
    InvalidArgument,
    RecordNameChangeIntegrityError,
    ValidationError,
)

if TYPE_CHECKING:
    import pandas as pd

    from lamindb.base.types import StrField


IPYTHON = getattr(builtins, "__IPYTHON__", False)


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


def validate_fields(record: Record, kwargs):
    from lamindb.base.validation import validate_literal_fields

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

        logger.warning(f"{msg}")
        if settings._verbosity_int >= 1:
            display(queryset.df())
    else:
        logger.warning(f"{msg}\n{queryset}")
    return None


def __init__(record: Record, *args, **kwargs):
    skip_validation = kwargs.pop("_skip_validation", False)
    if not args and skip_validation:
        super(BasicRecord, record).__init__(**kwargs)
    elif not args and not skip_validation:
        validate_fields(record, kwargs)

        # do not search for names if an id is passed; this is important
        # e.g. when synching ids from the notebook store to lamindb
        has_consciously_provided_uid = False
        if "_has_consciously_provided_uid" in kwargs:
            has_consciously_provided_uid = kwargs.pop("_has_consciously_provided_uid")
        if (
            isinstance(record, (CanCurate, Collection, Transform))
            and settings.creation.search_names
            and not has_consciously_provided_uid
        ):
            name_field = getattr(record, "_name_field", "name")
            exact_match = suggest_records_with_similar_names(record, name_field, kwargs)
            if exact_match is not None:
                if "version" in kwargs:
                    if kwargs["version"] is not None:
                        version_comment = " and version"
                        existing_record = record.__class__.filter(
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
                        f"returning existing {record.__class__.__name__} record with same"
                        f" {name_field}{version_comment}: '{kwargs[name_field]}'"
                    )
                    if isinstance(record, Schema):
                        if existing_record.hash != kwargs["hash"]:
                            raise ValueError(
                                f"Schema name is already in use by schema with uid '{existing_record.uid}', please choose a different name."
                            )
                    init_self_from_db(record, existing_record)
                    update_attributes(record, kwargs)
                    return None
        super(BasicRecord, record).__init__(**kwargs)
        if isinstance(record, ValidateFields):
            # this will trigger validation against django validators
            try:
                if hasattr(record, "clean_fields"):
                    record.clean_fields()
                else:
                    record._Model__clean_fields()
            except DjangoValidationError as e:
                message = _format_django_validation_error(record, e)
                raise FieldValidationError(message) from e
    elif len(args) != len(record._meta.concrete_fields):
        raise FieldValidationError(
            f"Use keyword arguments instead of positional arguments, e.g.: {record.__class__.__name__}(name='...')."
        )
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(BasicRecord, record).__init__(*args, **kwargs)
        _store_record_old_name(record)
        _store_record_old_key(record)


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

    for overload in overloads:
        params_block = overload.group(1)
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


@classmethod  # type:ignore
@doc_args(Record.filter.__doc__)
def filter(cls, *queries, **expressions) -> QuerySet:
    """{}"""  # noqa: D415
    from lamindb._query_set import QuerySet

    _using_key = None
    if "_using_key" in expressions:
        _using_key = expressions.pop("_using_key")

    return QuerySet(model=cls, using=_using_key).filter(*queries, **expressions)


@classmethod  # type:ignore
@doc_args(Record.get.__doc__)
def get(
    cls,
    idlike: int | str | None = None,
    **expressions,
) -> Record:
    """{}"""  # noqa: D415
    from lamindb._query_set import QuerySet

    return QuerySet(model=cls).get(idlike, **expressions)


@classmethod  # type:ignore
@doc_args(Record.df.__doc__)
def df(
    cls,
    include: str | list[str] | None = None,
    features: bool | list[str] = False,
    limit: int = 100,
) -> pd.DataFrame:
    """{}"""  # noqa: D415
    query_set = cls.filter()
    if hasattr(cls, "updated_at"):
        query_set = query_set.order_by("-updated_at")
    return query_set[:limit].df(include=include, features=features)


def _search(
    cls,
    string: str,
    *,
    field: StrField | list[StrField] | None = None,
    limit: int | None = 20,
    case_sensitive: bool = False,
    using_key: str | None = None,
    truncate_string: bool = False,
) -> QuerySet:
    if string is None:
        raise ValueError("Cannot search for None value! Please pass a valid string.")

    input_queryset = _queryset(cls, using_key=using_key)
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


@classmethod  # type: ignore
@doc_args(Record.search.__doc__)
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


def _lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
    using_key: str | None = None,
) -> NamedTuple:
    """{}"""  # noqa: D415
    queryset = _queryset(cls, using_key=using_key)
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


@classmethod  # type: ignore
@doc_args(Record.lookup.__doc__)
def lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
) -> NamedTuple:
    """{}"""  # noqa: D415
    return _lookup(cls=cls, field=field, return_field=return_field)


def get_name_field(
    registry: type[Record] | QuerySet | Manager,
    *,
    field: str | StrField | None = None,
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


def _queryset(cls: Record | QuerySet | Manager, using_key: str) -> QuerySet:
    if isinstance(cls, (QuerySet, Manager)):
        return cls.all()
    elif using_key is None or using_key == "default":
        return cls.objects.all()
    else:
        # using must be called on cls, otherwise the connection isn't found
        return cls.using(using_key).all()


def add_db_connection(db: str, using: str):
    db_config = dj_database_url.config(
        default=db, conn_max_age=600, conn_health_checks=True
    )
    db_config["TIME_ZONE"] = "UTC"
    db_config["OPTIONS"] = {}
    db_config["AUTOCOMMIT"] = True
    connections.settings[using] = db_config


@classmethod  # type: ignore
@doc_args(Record.using.__doc__)
def using(
    cls,
    instance: str | None,
) -> QuerySet:
    """{}"""  # noqa: D415
    from ._query_set import QuerySet

    if instance is None:
        return QuerySet(model=cls, using=None)
    owner, name = get_owner_name_from_identifier(instance)
    settings_file = instance_settings_file(name, owner)
    cache_filepath = ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    if not settings_file.exists():
        result = connect_instance_hub(owner=owner, name=name)
        if isinstance(result, str):
            raise RuntimeError(
                f"Failed to load instance {instance}, please check your permissions!"
            )
        iresult, _ = result
        source_module = {
            modules for modules in iresult["schema_str"].split(",") if modules != ""
        }  # type: ignore
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
    from lamindb.core._context import context
    from lamindb.core._data import WARNING_RUN_TRANSFORM

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


# docstring handled through attach_func_to_class_method
def save(self, *args, **kwargs) -> Record:
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
    transfer_logs: dict[str, list[str]] = {"mapped": [], "transferred": [], "run": None}
    if db is not None and db != "default" and using_key is None:
        if isinstance(self, IsVersioned):
            if not self.is_latest:
                raise NotImplementedError(
                    "You are attempting to transfer a record that's not the latest in its version history. This is currently not supported."
                )
        pre_existing_record = transfer_to_default_db(
            self, using_key, transfer_logs=transfer_logs
        )
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
                    super(BasicRecord, self).save(*args, **kwargs)  # type: ignore
                self._revises = None
            # save unversioned record
            else:
                super(BasicRecord, self).save(*args, **kwargs)
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
        _store_record_old_name(self)
        _store_record_old_key(self)
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

            from lamindb.models import FeatureManager

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
    return self


def _store_record_old_name(record: Record):
    # writes the name to the _name attribute, so we can detect renaming upon save
    if hasattr(record, "_name_field"):
        record._old_name = getattr(record, record._name_field)


def _store_record_old_key(record: Record):
    # writes the key to the _old_key attribute, so we can detect key changes upon save
    if isinstance(record, (Artifact, Transform)):
        record._old_key = record.key


def check_name_change(record: Record):
    """Warns if a record's name has changed."""
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


def check_key_change(record: Artifact | Transform):
    """Errors if a record's key has falsely changed."""
    if not isinstance(record, Artifact) or not hasattr(record, "_old_key"):
        return

    old_key = record._old_key or ""
    new_key = record.key or ""

    if old_key != new_key:
        if not record._key_is_virtual:
            raise InvalidArgument(
                f"Changing a non-virtual key of an artifact is not allowed! Tried to change key from '{old_key}' to '{new_key}'."
            )
        old_key_suffix = (
            record.suffix
            if record.suffix
            else extract_suffix_from_path(PurePosixPath(old_key), arg_name="key")
        )
        new_key_suffix = extract_suffix_from_path(
            PurePosixPath(new_key), arg_name="key"
        )
        if old_key_suffix != new_key_suffix:
            raise InvalidArgument(
                f"The suffix '{new_key_suffix}' of the provided key is incorrect, it should be '{old_key_suffix}'."
            )


def delete(self) -> None:
    """Delete the record."""
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
                super(BasicRecord, self).delete()  # type: ignore
            logger.warning(f"new latest version is {new_latest}")
            return None
    super(BasicRecord, self).delete()


METHOD_NAMES = [
    "__init__",
    "filter",
    "get",
    "df",
    "search",
    "lookup",
    "save",
    "delete",
    "using",
]

if ln_setup._TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(Record, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, BasicRecord, globals())
    attach_func_to_class_method(name, Record, globals())
