import builtins
from typing import Dict, Iterable, List, Literal, NamedTuple, Optional, Set, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_utils import colors, logger
from lamin_utils._lookup import Lookup
from lamin_utils._search import search as base_search
from lamindb_setup.dev._docs import doc_args
from lnschema_core import ORM
from lnschema_core.models import format_datetime
from lnschema_core.types import ListLike, StrField

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._feature_manager import create_features_df
from ._from_values import _has_species_field, get_or_create_records
from .dev._settings import settings

IPYTHON = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def init_self_from_db(self: ORM, existing_record: ORM):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def validate_required_fields(orm: ORM, kwargs):
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


def suggest_objects_with_same_name(orm: ORM, kwargs) -> Optional[str]:
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
                msg = (
                    "Records with similar names exist! Did you mean to load one of"
                    " them?"
                )
                if IPYTHON:
                    from IPython.display import display

                    logger.warning(f"{msg}")
                    display(results)
                else:
                    logger.warning(f"{msg}\n{results}")
    return None


def __init__(orm: ORM, *args, **kwargs):
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
                        f"Loaded record with exact same name{version_comment}"
                    )
                    init_self_from_db(orm, existing_record)
                    return None
        super(ORM, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("Please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(ORM, orm).__init__(*args, **kwargs)


def view_parents(
    self,
    field: Optional[StrField] = None,
    with_children: bool = False,
    distance: int = 100,
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
@doc_args(ORM.from_values.__doc__)
def from_values(cls, values: ListLike, field: StrField, **kwargs) -> List["ORM"]:
    """{}"""
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an ORM field, e.g., `CellType.name`!"
        )
    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    return get_or_create_records(
        iterable=values, field=field, from_bionty=from_bionty, **kwargs
    )


# From: https://stackoverflow.com/a/37648265
def _order_queryset_by_ids(queryset: QuerySet, ids: Iterable):
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
        return _order_queryset_by_ids(query_set, result.reset_index()["id"])
    else:
        return result.fillna("")


@classmethod  # type: ignore
@doc_args(ORM.search.__doc__)
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
@doc_args(ORM.lookup.__doc__)
def lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    """{}"""
    return _lookup(cls=cls, field=field)


def _inspect(
    cls,
    identifiers: ListLike,
    field: StrField,
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    from lamin_utils._inspect import inspect

    if not isinstance(field, str):
        field = field.field.name

    cls = cls.model if isinstance(cls, QuerySet) else cls

    return inspect(
        df=_filter_df_based_on_species(orm=cls, species=kwargs.get("species")),
        identifiers=identifiers,
        field=str(field),
        case_sensitive=case_sensitive,
        inspect_synonyms=inspect_synonyms,
        return_df=return_df,
        logging=logging,
    )


@classmethod  # type: ignore
@doc_args(ORM.inspect.__doc__)
def inspect(
    cls,
    identifiers: ListLike,
    field: StrField,
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    return _inspect(
        cls=cls,
        identifiers=identifiers,
        field=field,
        case_sensitive=case_sensitive,
        inspect_synonyms=inspect_synonyms,
        return_df=return_df,
        logging=logging,
        **kwargs,
    )


def _map_synonyms(
    cls,
    synonyms: Iterable,
    *,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    field: Optional[str] = None,
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    from lamin_utils._map_synonyms import map_synonyms

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    cls = cls.model if isinstance(cls, QuerySet) else cls

    try:
        cls._meta.get_field(synonyms_field)
        df = _filter_df_based_on_species(orm=cls, species=kwargs.get("species"))
    except FieldDoesNotExist:
        df = pd.DataFrame()
    return map_synonyms(
        df=df,
        identifiers=synonyms,
        field=field,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
    )


@classmethod  # type: ignore
@doc_args(ORM.map_synonyms.__doc__)
def map_synonyms(
    cls,
    synonyms: Iterable,
    *,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    field: Optional[str] = None,
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    return _map_synonyms(
        cls=cls,
        synonyms=synonyms,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
        field=field,
        **kwargs,
    )


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
    emojis = {"storage": "💾", "created_by": "👤", "transform": "💫", "run": "🚗"}
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
        values = [format_datetime(i.run_at) for i in self.input_of.all()]
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
    feature_sets = self.feature_sets.exclude(ref_field__startswith="core.Feature")
    if feature_sets.exists():
        for feature_set in feature_sets.all():
            key_split = feature_set.ref_field.split(".")
            if len(key_split) != 3:
                logger.warning(
                    "You have a legacy entry in feature_set.field, should be format"
                    " 'bionty.Gene.symbol'"
                )
            orm_name_with_schema = f"{key_split[0]}.{key_split[1]}"
            field_name = key_split[2]
            related_name = feature_sets_related_models.get(orm_name_with_schema)
            values = (
                feature_set.__getattribute__(related_name).all()[:5].list(field_name)
            )
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
    feature_sets = self.feature_sets.filter(ref_field__startswith="core.Feature").all()
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
                labels = self.features.get_labels(row["name"], mute=True)
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
    settings.verbosity = 2
    logger.info(msg)
    settings.verbosity = verbosity


def set_abbr(self, value: str):
    try:
        self.add_synonym(value, save=False)
    except NotImplementedError:
        pass
    self.abbr = value
    if not self._state.adding:
        self.save()


def _filter_df_based_on_species(
    orm: Union[ORM, QuerySet], species: Optional[Union[str, ORM]] = None
):
    import pandas as pd

    records = orm.all() if isinstance(orm, QuerySet) else orm.objects.all()
    if _has_species_field(orm):
        # here, we can safely import lnschema_bionty
        from lnschema_bionty._bionty import create_or_get_species_record

        species_record = create_or_get_species_record(species=species, orm=orm)
        if species_record is not None:
            records = records.filter(species__name=species_record.name)

    return pd.DataFrame.from_records(records.values())


def get_default_str_field(orm: Union[ORM, QuerySet, Manager]) -> str:
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


def _add_or_remove_synonyms(
    synonym: Union[str, Iterable],
    record: ORM,
    action: Literal["add", "remove"],
    force: bool = False,
    save: Optional[bool] = None,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: Set[str], record: ORM):
        """Errors if input synonym is associated with other records in the DB."""
        import pandas as pd
        from IPython.display import display

        syns_all = (
            record.__class__.objects.exclude(synonyms="").exclude(synonyms=None).all()
        )
        if len(syns_all) == 0:
            return
        df = pd.DataFrame(syns_all.values())
        df["synonyms"] = df["synonyms"].str.split("|")
        df = df.explode("synonyms")
        matches_df = df[(df["synonyms"].isin(synonyms)) & (df["id"] != record.id)]
        if matches_df.shape[0] > 0:
            records_df = pd.DataFrame(syns_all.filter(id__in=matches_df["id"]).values())
            logger.error(
                f"Input synonyms {matches_df['synonyms'].unique()} already associated"
                " with the following records:\n(Pass `force=True` to ignore this error)"
            )
            display(records_df)
            raise SystemExit(AssertionError)

    # passed synonyms
    if isinstance(synonym, str):
        syn_new_set = set([synonym])
    else:
        syn_new_set = set(synonym)
    # nothing happens when passing an empty string or list
    if len(syn_new_set) == 0:
        return
    # because we use | as the separator
    if any(["|" in i for i in syn_new_set]):
        raise AssertionError("A synonym can't contain '|'!")

    # existing synonyms
    syns_exist = record.synonyms
    if syns_exist is None or len(syns_exist) == 0:
        syns_exist_set = set()
    else:
        syns_exist_set = set(syns_exist.split("|"))

    if action == "add":
        if not force:
            check_synonyms_in_all_records(syn_new_set, record)
        syns_exist_set.update(syn_new_set)
    elif action == "remove":
        syns_exist_set = syns_exist_set.difference(syn_new_set)

    if len(syns_exist_set) == 0:
        syns_str = None
    else:
        syns_str = "|".join(syns_exist_set)

    record.synonyms = syns_str

    if save is None:
        # if record is already in DB, save the changes to DB
        save = not record._state.adding
    if save:
        record.save()


def _check_synonyms_field_exist(record: ORM):
    try:
        record.__getattribute__("synonyms")
    except AttributeError:
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        )


def add_synonym(
    self,
    synonym: Union[str, ListLike],
    force: bool = False,
    save: Optional[bool] = None,
):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(
        synonym=synonym, record=self, force=force, action="add", save=save
    )


def remove_synonym(self, synonym: Union[str, ListLike]):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


METHOD_NAMES = [
    "__init__",
    "search",
    "lookup",
    "map_synonyms",
    "inspect",
    "add_synonym",
    "remove_synonym",
    "from_values",
    "describe",
    "set_abbr",
    "view_parents",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(ORM, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, ORM, globals())


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
    logger.warning("select() is deprecated! Please rename: ORM.filter()")
    return cls.filter(**expressions)


@classmethod  # type: ignore
def select(cls, **expressions):
    return select_backward(cls, **expressions)


setattr(ORM, "__get_schema_name__", __get_schema_name__)
setattr(ORM, "__get_name_with_schema__", __get_name_with_schema__)
setattr(ORM, "select", select)  # backward compat
