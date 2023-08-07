import builtins
from typing import Dict, Iterable, List, NamedTuple, Optional, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_utils import colors, logger
from lamin_utils._lookup import Lookup
from lamin_utils._search import search as base_search
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Registry
from lnschema_core.models import format_field_value
from lnschema_core.types import ListLike, StrField

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records
from .dev._feature_manager import create_features_df
from .dev._settings import settings
from .dev._view_parents import _transform_emoji

IPYTHON = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def init_self_from_db(self: Registry, existing_record: Registry):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def validate_required_fields(orm: Registry, kwargs):
    required_fields = {
        k.name for k in orm._meta.fields if not k.null and k.default is None
    }
    required_fields_not_passed = {k: None for k in required_fields if k not in kwargs}
    kwargs.update(required_fields_not_passed)
    missing_fields = [
        k for k, v in kwargs.items() if v is None and k in required_fields
    ]
    if missing_fields:
        raise TypeError(f"{missing_fields} are required.")


def suggest_objects_with_same_name(orm: Registry, kwargs) -> Optional[str]:
    if kwargs.get("name") is None:
        return None
    else:
        results = orm.search(kwargs["name"])
        if results.shape[0] == 0:
            return None

        # subset results to those with at least 0.85 levensteihn distance
        results = results.loc[results.__ratio__ >= 85]

        # test for exact match
        if len(results) > 0:
            if results.index[0] == kwargs["name"]:
                return "object-with-same-name-exists"
            else:
                s = "" if results.shape[0] == 1 else "s"
                it = "it" if results.shape[0] == 1 else "one of them"
                msg = (
                    f"record{s} with similar name{s} exist! did you mean to load {it}?"
                )
                if IPYTHON:
                    from IPython.display import display

                    logger.warning(f"{msg}")
                    display(results)
                else:
                    logger.warning(f"{msg}\n{results}")
    return None


def __init__(orm: Registry, *args, **kwargs):
    if not args:
        validate_required_fields(orm, kwargs)
        from .dev._settings import settings

        if settings.upon_create_search_names:
            result = suggest_objects_with_same_name(orm, kwargs)
            if result == "object-with-same-name-exists":
                if "version" in kwargs:
                    version_comment = " and version "
                    existing_record = orm.filter(
                        name=kwargs["name"], version=kwargs["version"]
                    ).one_or_none()
                else:
                    version_comment = " "
                    existing_record = orm.filter(name=kwargs["name"]).one()
                if existing_record is not None:
                    logger.success(
                        f"loaded record with exact same name{version_comment}"
                    )
                    init_self_from_db(orm, existing_record)
                    return None
        super(Registry, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(Registry, orm).__init__(*args, **kwargs)


def view_parents(
    self,
    field: Optional[StrField] = None,
    with_children: bool = False,
    distance: int = 5,
):
    from lamindb.dev._view_parents import view_parents as _view_parents

    if field is None:
        field = get_default_str_field(self)
    if not isinstance(field, str):
        field = field.field.name

    return _view_parents(
        record=self, field=field, with_children=with_children, distance=distance
    )


@classmethod  # type:ignore
@doc_args(Registry.from_values.__doc__)
def from_values(cls, values: ListLike, field: StrField, **kwargs) -> List["Registry"]:
    """{}"""
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an Registry field, e.g., `CellType.name`!"
        )
    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    return get_or_create_records(
        iterable=values, field=field, from_bionty=from_bionty, **kwargs
    )


# From: https://stackoverflow.com/a/37648265
def _order_query_set_by_ids(queryset: QuerySet, ids: Iterable):
    from django.db.models import Case, When

    preserved = Case(*[When(pk=pk, then=pos) for pos, pk in enumerate(ids)])
    return queryset.filter(pk__in=ids).order_by(preserved)


def _search(
    cls,
    string: str,
    *,
    field: Optional[Union[StrField, List[StrField]]] = None,
    return_queryset: bool = False,
    limit: Optional[int] = None,
    case_sensitive: bool = False,
    synonyms_field: Optional[StrField] = "synonyms",
) -> Union["pd.DataFrame", "QuerySet"]:
    query_set = cls.all() if isinstance(cls, QuerySet) else cls.objects.all()
    orm = cls.model if isinstance(cls, QuerySet) else cls

    def _search_single_field(
        string: str,
        field: Optional[StrField],
        synonyms_field: Optional[StrField] = "synonyms",
    ) -> "pd.DataFrame":
        if field is None:
            field = get_default_str_field(cls)
        if not isinstance(field, str):
            field = field.field.name

        try:
            orm._meta.get_field(synonyms_field)
            synonyms_field_exists = True
        except FieldDoesNotExist:
            synonyms_field_exists = False

        if synonyms_field is not None and synonyms_field_exists:
            df = pd.DataFrame(query_set.values("id", field, synonyms_field))
        else:
            df = pd.DataFrame(query_set.values("id", field))

        return base_search(
            df=df,
            string=string,
            field=field,
            limit=limit,
            synonyms_field=str(synonyms_field),
            case_sensitive=case_sensitive,
        )

    # search in both key and description fields for file
    if orm.__name__ == "File" and field is None:
        field = ["key", "description"]

    if not isinstance(field, List):
        field = [field]

    results = []
    for fd in field:
        result_field = _search_single_field(
            string=string, field=fd, synonyms_field=synonyms_field
        )
        results.append(result_field)
        # turn off synonyms search after the 1st field
        synonyms_field = None

    if len(results) > 1:
        result = (
            pd.concat([r.reset_index() for r in results], join="outer")
            .drop(columns=["index"], errors="ignore")
            .set_index("id")
        )
    else:
        result = results[0]

    # remove results that have __ratio__ 0
    if "__ratio__" in result.columns:
        result = result[result["__ratio__"] > 0].sort_values(
            "__ratio__", ascending=False
        )
        # move the __ratio__ to be the last column
        result["__ratio__"] = result.pop("__ratio__")

    if return_queryset:
        return _order_query_set_by_ids(query_set, result.reset_index()["id"])
    else:
        return result.fillna("")


@classmethod  # type: ignore
@doc_args(Registry.search.__doc__)
def search(
    cls,
    string: str,
    *,
    field: Optional[StrField] = None,
    return_queryset: bool = False,
    limit: Optional[int] = None,
    case_sensitive: bool = False,
    synonyms_field: Optional[StrField] = "synonyms",
) -> Union["pd.DataFrame", "QuerySet"]:
    """{}"""
    return _search(
        cls=cls,
        string=string,
        field=field,
        return_queryset=return_queryset,
        limit=limit,
        case_sensitive=case_sensitive,
        synonyms_field=synonyms_field,
    )


def _lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    """{}"""
    if field is None:
        if cls._meta.model.__name__ == "User":
            field = cls._meta.get_field("handle").name
        else:
            field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.all() if isinstance(cls, QuerySet) else cls.objects.all()
    cls = cls.model if isinstance(cls, QuerySet) else cls

    return Lookup(
        records=records,
        values=[i.get(field) for i in records.values()],
        tuple_name=cls.__name__,
        prefix="ln",
    ).lookup()


@classmethod  # type: ignore
@doc_args(Registry.lookup.__doc__)
def lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    """{}"""
    return _lookup(cls=cls, field=field)


def describe(self):
    model_name = colors.green(self.__class__.__name__)
    msg = ""

    def dict_related_model_to_related_name(orm):
        d: Dict = {
            i.related_model.__get_name_with_schema__(): i.related_name
            for i in orm._meta.related_objects
            if i.related_name is not None
        }
        d.update(
            {
                i.related_model.__get_name_with_schema__(): i.name
                for i in orm._meta.many_to_many
                if i.name is not None
            }
        )

        return d

    # Display the file record
    fields = self._meta.fields
    direct_fields = []
    foreign_key_fields = []
    for f in fields:
        if f.is_relation:
            foreign_key_fields.append(f.name)
        else:
            direct_fields.append(f.name)

    # Display Provenance
    # display line by line the foreign key fields
    emojis = {
        "storage": "🗃️",
        "created_by": "👤",
        "transform": _transform_emoji(self.transform),
        "run": "🚗",
    }
    if len(foreign_key_fields) > 0:
        record_msg = f"{model_name}({''.join([f'{i}={self.__getattribute__(i)}, ' for i in direct_fields])})"  # noqa
        msg += f"{record_msg.rstrip(', )')})\n\n"

        msg += f"{colors.green('Provenance')}:\n    "
        related_msg = "".join(
            [
                f"{emojis.get(i, '📎')} {i}: {self.__getattribute__(i)}\n    "
                for i in foreign_key_fields
            ]
        )
        msg += related_msg
    # input of
    if self.input_of.exists():
        values = [format_field_value(i.run_at) for i in self.input_of.all()]
        msg += f"⬇️ input_of ({colors.italic('core.Run')}): {values}\n    "
    msg = msg.rstrip("    ")

    if not self.feature_sets.exists():
        print(msg)
        return
    else:
        feature_sets_related_models = dict_related_model_to_related_name(
            self.feature_sets.first()
        )
    # Display Features by slot
    msg += f"{colors.green('Features')}:\n"
    # var
    feature_sets = self.feature_sets.exclude(registry="core.Feature")
    if feature_sets.exists():
        for feature_set in feature_sets.all():
            key_split = feature_set.registry.split(".")
            if len(key_split) == 3:
                logger.warning(
                    "you have a legacy entry in feature_set.field, should be just"
                    " 'bionty.Gene'"
                )
            orm_name_with_schema = f"{key_split[0]}.{key_split[1]}"
            field_name = "id"
            related_name = feature_sets_related_models.get(orm_name_with_schema)
            values = feature_set.__getattribute__(related_name).all()[:5].list("id")
            slots = self.feature_sets.through.objects.filter(
                file=self, feature_set=feature_set
            ).list("slot")
            for slot in slots:
                if slot == "var":
                    slot += " (X)"
                msg += f"  🗺️ {colors.bold(slot)}:\n"
                ref = colors.italic(f"{orm_name_with_schema}.{field_name}")
                msg += f"    🔗 index ({feature_set.n}, {ref}): {values}\n".replace(
                    "]", "...]"
                )

    # obs
    # Feature, combine all features into one dataframe
    feature_sets = self.feature_sets.filter(registry="core.Feature").all()
    if feature_sets.exists():
        features_df = create_features_df(
            file=self, feature_sets=feature_sets.all(), exclude=True
        )
        for slot in features_df["slot"].unique():
            df_slot = features_df[features_df.slot == slot]
            if slot == "obs":
                slot += " (metadata)"
            msg += f"  🗺️ {colors.bold(slot)}:\n"
            for _, row in df_slot.iterrows():
                labels = self.get_labels(row["name"], mute=True)
                indent = ""
                if isinstance(labels, dict):
                    msg += f"    🔗 {row['name']} ({row.registries})\n"
                    indent = "    "
                else:
                    labels = {row["registries"]: labels}
                for registry, labels in labels.items():
                    count_str = f"{len(labels)}, {colors.italic(f'{registry}')}"
                    try:
                        field = get_default_str_field(labels)
                    except ValueError:
                        field = "id"
                    values = labels.list(field)[:5]
                    msg_objects = (
                        f"{indent}    🔗 {row['name']} ({count_str}): {values}\n"
                    )
                    msg += msg_objects
    msg = msg.rstrip("\n")
    msg = msg.rstrip("Features:")
    verbosity = settings.verbosity
    settings.verbosity = 3
    logger.info(msg)
    settings.verbosity = verbosity


def get_default_str_field(orm: Union[Registry, QuerySet, Manager]) -> str:
    """Get the 1st char or text field from the orm."""
    if isinstance(orm, (QuerySet, Manager)):
        orm = orm.model
    model_field_names = [i.name for i in orm._meta.fields]

    # set default field
    if orm._meta.model.__name__ == "Run":
        field = orm._meta.get_field("created_at")
    elif "name" in model_field_names:
        # by default use the name field
        field = orm._meta.get_field("name")
    else:
        # first char or text field that doesn't contain "id"
        for i in orm._meta.fields:
            if "id" in i.name:
                continue
            if i.get_internal_type() in {"CharField", "TextField"}:
                field = i
                break

    # no default field can be found
    if field is None:
        raise ValueError("Please specify a field to search against!")

    return field.name


METHOD_NAMES = [
    "__init__",
    "search",
    "lookup",
    "from_values",
    "describe",
    "view_parents",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(Registry, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Registry, globals())


@classmethod  # type: ignore
def __get_schema_name__(cls) -> str:
    schema_module_name = cls.__module__.split(".")[0]
    schema_name = schema_module_name.replace("lnschema_", "")
    return schema_name


@classmethod  # type: ignore
def __get_name_with_schema__(cls) -> str:
    schema_name = cls.__get_schema_name__()
    return f"{schema_name}.{cls.__name__}"


def select_backward(cls, **expressions):
    logger.warning("select() is deprecated! please use: Registry.filter()")
    return cls.filter(**expressions)


@classmethod  # type: ignore
def select(cls, **expressions):
    return select_backward(cls, **expressions)


setattr(Registry, "__get_schema_name__", __get_schema_name__)
setattr(Registry, "__get_name_with_schema__", __get_name_with_schema__)
setattr(Registry, "select", select)  # backward compat
