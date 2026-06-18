from __future__ import annotations

import os
from collections.abc import Iterable
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Literal, overload

import pgtrigger
from django.conf import settings as django_settings
from django.db import models
from django.db.models import CASCADE, PROTECT
from lamin_utils import logger

from lamindb.base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    JSONField,
    TextField,
)
from lamindb.base.utils import class_and_instance_method, strict_classmethod
from lamindb.errors import FieldValidationError

from ..base.uids import base62_16
from .artifact import Artifact
from .can_curate import CanCurate
from .collection import Collection
from .feature import Feature, convert_to_pandas_dtype
from .has_parents import HasParents, _query_relatives
from .query_set import (
    QuerySet,
    encode_lamindb_fields_as_columns,
    get_default_branch_ids,
    reorder_subset_columns_in_df,
)
from .run import Run, TracksRun, TracksUpdates, User, current_run, current_user_id
from .sqlrecord import (
    BaseSQLRecord,
    Branch,
    HasType,
    IsLink,
    Space,
    SQLRecord,
    _get_record_kwargs,
    pop_space_branch_kwargs,
)
from .transform import Transform
from .ulabel import ULabel

if TYPE_CHECKING:
    import pandas as pd

    from ._feature_manager import FeatureManager
    from .block import RecordBlock
    from .project import Project, RecordProject, RecordReference, Reference
    from .query_manager import RelatedManager
    from .query_set import SQLRecordList
    from .schema import Schema


# keep docstring in sync with test_record_docstring_examples in test_record_basics.py
IMPORTS_UID = "W3WdiFRZTvTJajNp"
SCHEMA_IMPORTS_UID = "DGZkj4yhGWMJE5fu"


def get_type_schema_index(record_type: Record | None) -> Feature | None:
    """Return the index feature for a record type sheet, if configured."""
    if record_type is None or not record_type.is_type:
        return None
    schema = record_type.schema
    if schema is None:
        return None
    return schema.index


def is_schema_index_feature(schema: Schema | None, feature: Feature) -> bool:
    if schema is None:
        return False
    index_feature = schema.index
    return index_feature is not None and feature.uid == index_feature.uid


def validate_record_sheet_index_feature(index_feature: Feature) -> None:
    """Ensure a record-sheet index feature can be stored on `Record.name`."""
    if index_feature.dtype_as_str != "str":
        raise ValueError(
            f"schema index feature '{index_feature.name}' must have dtype str "
            f"because it is stored on Record.name, not {index_feature.dtype_as_str!r}"
        )


def validate_record_type_schema_index(schema: Schema | None) -> None:
    if schema is None:
        return
    index_feature = schema.index
    if index_feature is not None:
        validate_record_sheet_index_feature(index_feature)


def coerce_index_value_to_record_name(value: Any, feature: Feature) -> str | None:
    """Convert an index feature value to a `Record.name` string."""
    import pandas as pd

    if value is None or (pd.api.types.is_scalar(value) and pd.isna(value)):
        return None
    if isinstance(value, str):
        return value
    raise TypeError(
        f"index feature '{feature.name}' value must be a string, not {type(value).__name__}"
    )


def index_value_from_record_name(name: str | None, feature: Feature) -> str | None:
    """Convert `Record.name` back to the index feature value."""
    return name


def apply_index_feature_to_record(
    record: Record,
    feature: Feature,
    value: Any,
    *,
    persist: bool = True,
) -> None:
    """Set `record.name` from an index feature value."""
    record.name = coerce_index_value_to_record_name(value, feature)
    if persist and record.pk is not None:
        persist_record_name(record)


def persist_record_name(record: Record) -> None:
    """Persist `Record.name` without re-entering lazy `features` saving."""
    SQLRecord.save(record, update_fields=["name"])


def inject_index_into_feature_dict(record: Record, dictionary: dict[str, Any]) -> None:
    """Expose the index feature in `get_values` / feature dicts from `record.name`."""
    index_feature = get_type_schema_index(record.type)
    if index_feature is None or record.name is None:
        return
    dictionary[index_feature.name] = index_value_from_record_name(
        record.name, index_feature
    )


def pop_index_from_feature_dictionary(
    dictionary: dict[str, Any], schema: Schema
) -> tuple[str | None, dict[str, Any]]:
    """Extract index value for `record.name` and remove it from feature payload."""
    index_feature = schema.index
    if index_feature is None:
        return None, dictionary
    index_name = index_feature.name
    if index_name not in dictionary:
        return None, dictionary
    value = dictionary.pop(index_name)
    return coerce_index_value_to_record_name(value, index_feature), dictionary


def export_includes_record_metadata(schema: Schema | None) -> bool:
    """Whether sheet export includes encoded ``__lamindb_record_*`` columns."""
    return schema is None or schema.index is None


def drop_record_metadata_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Drop encoded Record metadata columns from an export dataframe."""
    metadata_cols = [
        col
        for col in df.columns
        if col.startswith("__lamindb_record_") and col.endswith("__")
    ]
    if metadata_cols:
        df = df.drop(columns=metadata_cols)
    encoded_id = encode_lamindb_fields_as_columns(Record, "id")
    if isinstance(encoded_id, str) and df.index.name == encoded_id:
        df = df.reset_index(drop=True)
    return df


def apply_schema_index_to_export_dataframe(
    df: pd.DataFrame,
    index_feature: Feature,
    *,
    encoded_id: str,
    encoded_name: str,
    include_record_metadata: bool = True,
) -> pd.DataFrame:
    """Move the schema index feature from columns to `DataFrame.index`."""
    index_col = index_feature.name
    lamin_record_ids = df.index.to_series()
    if index_col in df.columns:
        index_values = df[index_col]
        df = df.drop(columns=[index_col])
    elif "name" in df.columns:
        index_values = df["name"]
        df = df.drop(columns=["name"])
    elif encoded_name in df.columns:
        index_values = df[encoded_name]
        df = df.drop(columns=[encoded_name])
    else:
        raise ValueError(
            f"could not find values for schema index feature '{index_col}' in export dataframe"
        )

    df = df.copy()
    if include_record_metadata:
        df[encoded_id] = lamin_record_ids.values
    df = df.set_index(index_values)
    df.index.name = index_col
    if not include_record_metadata:
        df = drop_record_metadata_columns(df)
    return df


def dataframe_for_record_batch(
    df: pd.DataFrame, index_feature: Feature | None
) -> pd.DataFrame:
    """Normalize a batch dataframe so the schema index is a column for row iteration."""
    import pandas as pd

    if index_feature is None:
        return df
    if df.index.name == index_feature.name and not isinstance(df.index, pd.RangeIndex):
        return df.reset_index()
    return df


def move_schema_index_column_to_dataframe_index(
    df: pd.DataFrame, schema: Schema
) -> pd.DataFrame:
    """Align a validation dataframe with `Schema.index` (index on `df.index`)."""
    index_feature = schema.index
    if index_feature is None or index_feature.name not in df.columns:
        return df
    df = df.copy()
    df = df.set_index(index_feature.name)
    df.index.name = index_feature.name
    return df


def strip_index_for_record_persistence(
    record: Record,
    schema: Schema,
    dictionary: dict[str, Any],
    feature_objects: list[Feature],
    *,
    values_by_feature_uid: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], list[Feature]]:
    """Move schema index values to `record.name` and drop them from link-table writes."""
    index_feature = schema.index
    if index_feature is None:
        return dictionary, feature_objects

    index_value = None
    if values_by_feature_uid is not None and index_feature.uid in values_by_feature_uid:
        index_value = values_by_feature_uid[index_feature.uid]
    elif index_feature.name in dictionary:
        index_value = dictionary[index_feature.name]

    if index_value is not None:
        apply_index_feature_to_record(record, index_feature, index_value, persist=False)

    dictionary = dict(dictionary)
    dictionary.pop(index_feature.name, None)
    feature_objects = [
        feature for feature in feature_objects if feature.uid != index_feature.uid
    ]
    return dictionary, feature_objects


IndexNameConflict = Literal["keep_name", "keep_feature"]


def _record_sheet_label(record: Record) -> str:
    if record.type_id is None:
        return "unknown sheet"
    sheet = record.type
    if sheet is None:
        return f"type_id={record.type_id}"
    if sheet.name:
        return f"{sheet.name} ({sheet.uid})"
    return sheet.uid


def _collect_promote_index_name_conflicts(
    json_links: list[tuple[int, int, Any]],
    records_by_id: dict[int, Record],
    feature: Feature,
) -> list[tuple[int, str, str, str]]:
    conflicts: list[tuple[int, str, str, str]] = []
    for _, record_id, value in json_links:
        record = records_by_id.get(record_id)
        if record is None:
            continue
        feature_name = coerce_index_value_to_record_name(value, feature)
        if feature_name is not None and record.name and record.name != feature_name:
            conflicts.append(
                (record_id, _record_sheet_label(record), record.name, feature_name)
            )
    return conflicts


def _resolve_index_name_conflict(
    conflicts: list[tuple[int, str, str, str]],
    feature: Feature,
    resolution: IndexNameConflict | None,
) -> IndexNameConflict:
    if not conflicts:
        return "keep_name"
    if resolution is not None:
        return resolution
    if os.getenv("LAMIN_TESTING") == "true":
        return "keep_name"

    feature_name = feature.name
    n = len(conflicts)
    example_lines = "\n".join(
        f"  sheet {sheet_label}, record {record_id}: Record.name={name!r}, {feature_name}={feature_value!r}"
        for record_id, sheet_label, name, feature_value in conflicts[:3]
    )
    if n > 3:
        example_lines += f"\n  ... and {n - 3} more"
    response = input(
        f"{n} sheet row(s) have both Record.name and '{feature_name}' values that differ.\n"
        f"{example_lines}\n"
        "Keep Record.name (y) or use feature values (n)? "
    )
    if response == "y":
        return "keep_name"
    if response == "n":
        return "keep_feature"
    raise ValueError("schema save cancelled: index name conflict not resolved")


def migrate_record_sheet_index_on_schema_save(
    schema: Schema,
    *,
    old_index_uid: str | None,
    new_index_uid: str | None,
    using: str | None = None,
    index_name_conflict: IndexNameConflict | None = None,
) -> None:
    """Migrate sheet row keys when ``schema.index`` changes on save."""
    if (
        old_index_uid == new_index_uid
        or schema.is_type
        or schema.itype != "Feature"
        or schema.otype == "Composite"
    ):
        return

    from django.db import transaction

    from .save import bulk_create

    db = using or schema._state.db
    sheet_type_ids = (
        Record.objects.using(db)
        .filter(is_type=True, schema_id=schema.id)
        .values_list("id", flat=True)
    )
    if not sheet_type_ids:
        return

    row_ids = (
        Record.objects.using(db)
        .filter(is_type=False, type_id__in=sheet_type_ids)
        .values_list("id", flat=True)
    )

    with transaction.atomic(using=db):
        if old_index_uid is not None:
            old_feature = Feature.objects.using(db).get(uid=old_index_uid)
            rows_with_name = (
                Record.objects.using(db)
                .filter(
                    id__in=row_ids,
                    name__isnull=False,
                )
                .exclude(name="")
            )
            json_to_create = [
                RecordJson(record_id=record_id, feature=old_feature, value=name)
                for record_id, name in rows_with_name.values_list("id", "name")
                if name
            ]
            if json_to_create:
                bulk_create(json_to_create, ignore_conflicts=True)
            rows_with_name.update(name=None)

        if new_index_uid is not None:
            new_feature = Feature.objects.using(db).get(uid=new_index_uid)
            validate_record_sheet_index_feature(new_feature)
            json_links = (
                RecordJson.objects.using(db)
                .filter(
                    record_id__in=row_ids,
                    feature=new_feature,
                )
                .values_list("id", "record_id", "value")
            )
            if not json_links:
                return

            records_by_id = {
                record.id: record
                for record in Record.objects.using(db)
                .filter(id__in={record_id for _, record_id, _ in json_links})
                .select_related("type")
            }
            conflicts = _collect_promote_index_name_conflicts(
                list(json_links), records_by_id, new_feature
            )
            resolution = _resolve_index_name_conflict(
                conflicts, new_feature, index_name_conflict
            )
            records_to_update: list[Record] = []
            json_ids_to_delete: list[int] = []
            for link_id, record_id, value in json_links:
                record = records_by_id.get(record_id)
                if record is None:
                    json_ids_to_delete.append(link_id)
                    continue
                name = coerce_index_value_to_record_name(value, new_feature)
                if name is not None:
                    if not record.name:
                        record.name = name
                        records_to_update.append(record)
                    elif record.name != name and resolution == "keep_feature":
                        record.name = name
                        records_to_update.append(record)
                json_ids_to_delete.append(link_id)

            if records_to_update:
                Record.objects.using(db).bulk_update(records_to_update, ["name"])
            if json_ids_to_delete:
                RecordJson.objects.using(db).filter(id__in=json_ids_to_delete).delete()


class RecordBatch:
    """DataFrame-backed batch created by :meth:`Record.from_dataframe`."""

    def __init__(
        self,
        *,
        cls: type[Record],
        df: pd.DataFrame,
        resolved_type: Record,
        name_field: str,
    ) -> None:
        self._cls = cls
        self._df = df
        self._resolved_type = resolved_type
        self._name_field = name_field
        self._records: list[Record] | None = None

    def __len__(self) -> int:
        return len(self._df)

    @property
    def type(self) -> Record:
        return self._resolved_type

    def _build_records(self) -> list[Record]:
        import pandas as pd

        index_feature = get_type_schema_index(self._resolved_type)
        records: list[Record] = []
        work_df = dataframe_for_record_batch(self._df, index_feature)
        row_dicts = work_df.to_dict(orient="records")
        for row in row_dicts:
            row = dict(row)
            name = None
            if index_feature is not None and index_feature.name in row:
                value = row.pop(index_feature.name)
                if not (pd.api.types.is_scalar(value) and pd.isna(value)):
                    name = coerce_index_value_to_record_name(value, index_feature)
            if name is None:
                if self._name_field in row:
                    name = row.pop(self._name_field)
                elif "name" in row:
                    name = row.pop("name")
                else:
                    name = None
                if pd.api.types.is_scalar(name) and pd.isna(name):
                    name = None

            features: dict[str, Any] = {}
            for key, value in row.items():
                if pd.api.types.is_scalar(value) and pd.isna(value):
                    continue
                features[key] = value

            if index_feature is not None and self._resolved_type.schema is not None:
                name_from_features, features = pop_index_from_feature_dictionary(
                    features, self._resolved_type.schema
                )
                if name is None:
                    name = name_from_features

            record_kwargs: dict[str, Any] = {"type": self._resolved_type}
            if features:
                record_kwargs["features"] = features
            records.append(self._cls(name=name, **record_kwargs))
        return records

    def save(self) -> SQLRecordList[Record]:
        """Persist all records and their feature values."""
        from .query_set import SQLRecordList
        from .save import save as ln_save

        if self._records is None:
            self._records = self._build_records()
        ln_save(self._records)
        return SQLRecordList(self._records)


class Record(SQLRecord, HasType, HasParents, CanCurate, TracksRun, TracksUpdates):
    """Flexible records with markdown notes, dynamic registries, and sheets.

    Useful for managing notes, experiments, samples, donors, cells, compounds, sequences,
    and other custom entities with their features.

    Args:
        name: `str | None = None` A name.
        description: `str | None = None` A description.
        type: `Record | None = None` The type of this record.
        is_type: `bool = False` Whether this record is a type (a record that
            classifies other records).
        features: `dict[str | Feature, Any] | None = None` Lazy feature values
            to persist on `.save()` or `ln.save([...])`.
        schema: `Schema | None = None` A schema defining allowed features for records of this type. Only applicable when `is_type=True`.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.
        branch: `Branch | None = None` A branch. If None, uses the current branch.
        space: `Space | None = None` A space. If None, uses the current space.


    See Also:
        :class:`~lamindb.Feature`
            Measurable properties such as columns of a sheet.
        :class:`~lamindb.Schema`
            Constrain sheet columns; :attr:`~lamindb.Schema.index` defines row keys.
        :class:`~lamindb.ULabel`
            Simple universal labels.

    Examples
    --------

    Create a **record** with a single feature::

        # create a feature if you don't yet have one
        gc_content = ln.Feature(name="gc_content", dtype=float).save()

        # create a record to track a sample
        sample1 = ln.Record(name="Sample 1", features={"gc_content": 0.5}).save()

        # describe the record
        sample1.describe()

    Group records in a dynamic registry by creating a **record type**, optionally constrained with a :class:`~lamindb.Schema`::

        # create an experiments registry
        experiments_registry = ln.Record(name="Experiments", is_type=True).save()
        experiment1 = ln.Record(name="Experiment 1", type=experiments_registry).save()

        # create a feature to link experiments
        experiment = ln.Feature(name="experiment", dtype=experiments_registry).save()

        # constrain a samples registry with a schema, turning it into a sheet
        schema = ln.Schema([experiment, gc_content.with_config(optional=True)], name="sample_schema").save()
        sample_sheet = ln.Record(name="Sample Sheet", is_type=True, schema=schema).save()

        # move the sample1 record into the sample sheet
        sample1.type = sample_sheet
        sample1.save()

        # reset the feature values for the record including the experiment
        sample1.features.set_values({
            gc_content: 0.5,
            experiment: "Experiment 1",  # automatically resolves by name, also accepts the experiment1 object
        })

    Export all records of a type to a dataframe::

        experiments_registry.to_dataframe()
        #> __lamindb_record_name__   ...
        #>            Experiment 1   ...
        #>            Experiment 2   ...

    Use :attr:`~lamindb.Schema.index` on a sheet schema to define row keys::

        sample_id = ln.Feature(name="sample_id", dtype=str).save()
        score = ln.Feature(name="score", dtype=float).save()
        schema = ln.Schema(features=[score], index=sample_id).save()
        sheet = ln.Record(name="Samples", is_type=True, schema=schema).save()

        record = ln.Record(type=sheet, features={"sample_id": "S-001", "score": 1.5}).save()
        assert record.name == "S-001"

        df = sheet.to_dataframe()
        assert df.index.name == "sample_id"
        assert "sample_id" not in df.columns

    Import records from a dataframe :meth:`~lamindb.Record.from_dataframe`::

        records = ln.Record.from_dataframe(df, type="my_df").save()  # creates a type my_df with inferred schema

    If you try to set incomplete features in a record in a sheet, you'll get a validation error::

        sample2 = ln.Record(name="Sample 2", type=sample_sheet).save()
        sample2.features.set_values({gc_content: 0.6})  # raises ValidationError because experiment is missing

    Query records by features:

    .. tab-set::

        .. tab-item:: Via objects

            .. code-block:: python

                ln.Record.filter(gc_content == 0.55)  # exact match
                ln.Record.filter(gc_content > 0.5)    # greater than

        .. tab-item:: Via strings

            .. code-block:: python

                ln.Record.filter(gc_content=0.55)  # exact match
                ln.Record.filter(gc_content__gt=0.5)    # greater than

    Query records by field::

        ln.Record.filter(type=sample_sheet)   # just the record on the sheet

    Notes
    -----

    **Schema index.** When a sheet schema defines :attr:`~lamindb.Schema.index`, the
    index feature acts as the row key — analogous to `df.index` for tabular data.
    The index feature must have `dtype=str` because values are stored on
    :attr:`~lamindb.Record.name`:

    - **Write**: `Record(features=...)`, `features.add_values()`, and
      :meth:`~lamindb.Record.from_dataframe` route the index feature to
      :attr:`~lamindb.Record.name` and do not write it to link tables.
    - **Read**: `features.get_values()` injects the index from `Record.name`.
    - **Export**: :meth:`~lamindb.Record.to_dataframe` puts the index on `df.index`
      (named after the index feature) and omits encoded metadata columns
      (`__lamindb_record_id__`, `__lamindb_record_uid__`, `__lamindb_record_name__`, etc.).
      Sheets without `index` keep the previous export behavior.
    - **Import**: :meth:`~lamindb.Record.from_dataframe` accepts a dataframe whose index
      matches the schema index feature (or the index feature as a column).
    - **CSV**: :meth:`~lamindb.Record.to_artifact` writes with `index=True` when an index
      is configured.

    You can edit records like spreadsheets on the hub:

    .. image:: https://lamin-site-assets.s3.amazonaws.com/.lamindb/XSzhWUb0EoHOejiw0002.png
        :width: 800px

    Just like for :class:`~lamindb.ULabel`, you can also model **ontologies** through the `parents`/`children` attributes.

    .. dropdown:: What is the difference between `Record` and `SQLRecord`?

        The features of a `Record` are flexible: you can dynamically define features and add features to a record.
        The fields of a `SQLRecord` are static: you need to define them in code and then migrate the underlying database.

        In complete analogy to this: A **record type** can model a registry dynamically, whereas a :class:`~lamindb.models.Registry` has to be
        defined as a static Python class together with its SQL database migration:  `lamin migrate create` and `lamin migrate deploy`.

        See :class:`~lamindb.models.SQLRecord` or the glossary for more information: :term:`docs:record`.

    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False
        app_label = "lamindb"
        if (
            django_settings.DATABASES.get("default", {}).get("ENGINE")
            == "django.db.backends.postgresql"
        ):
            triggers = [
                pgtrigger.Trigger(
                    name="prevent_record_type_cycle",
                    operation=pgtrigger.Update | pgtrigger.Insert,
                    when=pgtrigger.Before,
                    condition=pgtrigger.Condition("NEW.type_id IS NOT NULL"),
                    func="""
                        -- Enforce that records typed by a single-space type stay in the same space
                        IF EXISTS (
                            SELECT 1
                            FROM lamindb_record r
                            WHERE r.id = NEW.type_id
                              AND (
                                (
                                  r._aux->>'ss' = '1'
                                  AND r.space_id IS DISTINCT FROM NEW.space_id
                                )
                                OR (
                                  r._aux ? 'ss'
                                  AND r._aux->>'ss' <> '1'
                                  AND (
                                      SELECT sp.uid
                                      FROM lamindb_space sp
                                      WHERE sp.id = NEW.space_id
                                  ) IS DISTINCT FROM r._aux->>'ss'
                                )
                              )
                        ) THEN
                            RAISE EXCEPTION 'Cannot set type: record space must match locked type space';
                        END IF;

                        -- Check for direct self-reference
                        IF NEW.type_id = NEW.id THEN
                            RAISE EXCEPTION 'Cannot set type: record cannot be its own type';
                        END IF;

                        -- Check for cycles in the type chain
                        IF EXISTS (
                            WITH RECURSIVE type_chain AS (
                                SELECT type_id, 1 as depth
                                FROM lamindb_record
                                WHERE id = NEW.type_id

                                UNION ALL

                                SELECT r.type_id, tc.depth + 1
                                FROM lamindb_record r
                                INNER JOIN type_chain tc ON r.id = tc.type_id
                                WHERE tc.depth < 100
                            )
                            SELECT 1 FROM type_chain WHERE type_id = NEW.id
                        ) THEN
                            RAISE EXCEPTION 'Cannot set type: would create a cycle';
                        END IF;

                        RETURN NEW;
                    """,
                ),
            ]
        # also see raw SQL constraints for `is_type` and `type` FK validity in migrations

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, db_index=True, max_length=16, default=base62_16
    )
    """A universal random id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True, null=True)
    """Name or title of record (optional)."""
    type: Record | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of record, e.g., `Sample`, `Donor`, `Cell`, `Compound`, `Sequence` ← :attr:`~lamindb.Record.records`.

    Allows to group records by type, e.g., all samples, all donors, all cells, all compounds, all sequences.
    """
    records: RelatedManager[Record]
    """If a `type` (`is_type=True`), records of this `type`."""
    description: str | None = TextField(null=True)
    """A description."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A simple reference like a URL or external ID."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of simple reference."""
    extra_data: dict | None = models.JSONField(null=True)
    """Extra data in JSON format, not validated as features."""
    schema: Schema | None = ForeignKey(
        "Schema", CASCADE, null=True, related_name="records"
    )
    """A schema to enforce for a type ← :attr:`~lamindb.Schema.records`.

    This is analogous to the `schema` attribute of an `Artifact`.
    If `is_type` is `True`, the schema is used to enforce features for each record of this type.
    """
    linked_records: RelatedManager[Record] = models.ManyToManyField(
        "Record",
        through="RecordRecord",
        symmetrical=False,
        related_name="linked_in_records",
    )
    """Records linked in this record as a value ← :attr:`~lamindb.Record.linked_in_records`."""
    linked_in_records: RelatedManager[Record]
    """Records linking this record as a value. Is reverse accessor for `linked_records`."""
    parents: RelatedManager[Record] = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Ontological parents of this record ← :attr:`~lamindb.Record.children`.

    You can build an ontology under a given `type`. For example, introduce a type `CellType` and model the hiearchy of cell types under it via `parents` and `children`.
    """
    children: RelatedManager[Record]
    """Ontological children of this record. Is reverse accessor for `parents`."""
    # this is handled manually here because we want to se the related_name attribute
    # (this doesn't happen via inheritance of TracksRun, everything else is the same)
    run: Run | None = ForeignKey(
        Run,
        PROTECT,
        related_name="output_records",
        null=True,
        default=current_run,
        editable=False,
    )
    """Run that created the record ← :attr:`~lamindb.Run.output_records`."""
    input_of_runs: RelatedManager[Run] = models.ManyToManyField(
        Run, related_name="input_records"
    )
    """Runs that use this record as an input ← :attr:`~lamindb.Run.input_records`."""
    artifacts: RelatedManager[Artifact] = models.ManyToManyField(
        Artifact, through="ArtifactRecord", related_name="records"
    )
    """Artifacts annotated by this record ← :attr:`~lamindb.Artifact.records`."""
    runs: RelatedManager[Run] = models.ManyToManyField(
        Run, through="RunRecord", related_name="records"
    )
    """Runs annotated by this record ← :attr:`~lamindb.Run.records`."""
    transforms: RelatedManager[Transform] = models.ManyToManyField(
        Transform, through="TransformRecord", related_name="records"
    )
    """Transforms annotated by this record ← :attr:`~lamindb.Transform.records`."""
    collections: RelatedManager[Collection] = models.ManyToManyField(
        Collection, through="CollectionRecord", related_name="records"
    )
    """Collections annotated by this record ← :attr:`~lamindb.Collection.records`."""
    projects: RelatedManager[Project]
    """Projects that annotate this record ← :attr:`~lamindb.Project.records`."""
    references: RelatedManager[Reference]
    """References that annotate this record ← :attr:`~lamindb.Reference.records`."""
    linked_transforms: RelatedManager[Transform]
    """Transforms linked in this record as values ← :attr:`~lamindb.Transform.linked_in_records`."""
    linked_runs: RelatedManager[Run]
    """Runs linked in this record as values ← :attr:`~lamindb.Run.linked_in_records`."""
    linked_ulabels: RelatedManager[ULabel]
    """ULabels linked in this record as values ← :attr:`~lamindb.ULabel.linked_in_records`."""
    linked_artifacts: RelatedManager[Artifact]
    """Artifacts linked in this record as values ← :attr:`~lamindb.Artifact.linked_in_records`."""
    linked_projects: RelatedManager[Project]
    """Projects linked in this record as values ← :attr:`~lamindb.Project.linked_in_records`."""
    linked_references: RelatedManager[Reference]
    """References linked in this record as values ← :attr:`~lamindb.Reference.linked_in_records`."""
    linked_collections: RelatedManager[Collection]
    """Collections linked in this record as values ← :attr:`~lamindb.Collection.linked_in_records`."""
    linked_users: RelatedManager[User]
    """Users linked in this record as values ← :attr:`~lamindb.User.linked_in_records`."""
    ablocks: RelatedManager[RecordBlock]
    """Attached blocks ← :attr:`~lamindb.RecordBlock.record`."""
    values_json: RelatedManager[RecordJson]
    """JSON values `(record_id, feature_id, value)`."""
    values_record: RelatedManager[RecordRecord]
    """Record values with their features `(record_id, feature_id, value_id)`."""
    values_ulabel: RelatedManager[RecordULabel]
    """ULabel values with their features `(record_id, feature_id, value_id)`."""
    values_user: RelatedManager[RecordUser]
    """User values with their features `(record_id, feature_id, value_id)`."""
    values_transform: RelatedManager[RecordTransform]
    """Transform values with their features `(record_id, feature_id, value_id)`."""
    values_run: RelatedManager[RecordRun]
    """Run values with their features `(record_id, feature_id, value_id)`."""
    values_artifact: RelatedManager[RecordArtifact]
    """Artifact values with their features `(record_id, feature_id, value_id)`."""
    values_collection: RelatedManager[RecordCollection]
    """Collection values with their features `(record_id, feature_id, value_id)`."""
    values_reference: RelatedManager[RecordReference]
    """Reference values with their features `(record_id, feature_id, value_id)`."""
    values_project: RelatedManager[RecordProject]
    """Project values with their features `(record_id, feature_id, value_id)`."""

    @overload
    def __init__(
        self,
        name: str | None = None,
        type: Record | None = None,
        is_type: bool = False,
        features: dict[str | Feature, Any] | None = None,
        description: str | None = None,
        schema: Schema | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        branch: Branch | None = None,
        space: Space | None = None,
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
        if len(args) > 0:
            raise ValueError("Only one non-keyword arg allowed")
        name: str = kwargs.pop("name", None)
        type: str | None = kwargs.pop("type", None)
        is_type: bool = kwargs.pop("is_type", False)
        features: dict[str | Feature, Any] | None = kwargs.pop("features", None)
        description: str | None = kwargs.pop("description", None)
        schema: Schema | None = kwargs.pop("schema", None)
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        space_branch_kwargs = pop_space_branch_kwargs(kwargs)
        _skip_validation = kwargs.pop("_skip_validation", False)
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Record)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        if schema and not is_type:
            logger.important("passing schema, treating as type")
            is_type = True
        if features is not None:
            self._features = features
        super().__init__(
            name=name,
            type=type,
            is_type=is_type,
            description=description,
            reference=reference,
            reference_type=reference_type,
            schema=schema,
            _skip_validation=_skip_validation,
            _aux=_aux,
            **space_branch_kwargs,
        )

    def save(self, *args, **kwargs) -> Record:
        if self.is_type:
            validate_record_type_schema_index(self.schema)
        super().save(*args, **kwargs)
        if hasattr(self, "_features"):
            pending_features = self._features
            self.features.add_values(pending_features)
            del self._features
        return self

    @strict_classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        *,
        type: Record | str,
        name_field: str = "__lamindb_record_name__",
    ) -> RecordBatch:
        """Construct a dataframe-backed batch of records for bulk saving.

        Returns a :class:`RecordBatch`. Follow with `records.save()`.

        When the target sheet schema defines :attr:`~lamindb.Schema.index`, the index
        feature may be passed on `df.index` (named after the feature) or as a column.

        Args:
            df: A dataframe where rows represent records.
            type: Record type for all rows as either a `Record` object or a
                string. If passing a string, a new type with that name is created
                under `Imports` with an inferred schema from the dataframe.
                If that type name already exists, raise an error and pass an
                existing `Record` object for reuse.
                If the resolved type is a sheet (`type.schema is not None`), feature
                values are validated against that schema at save time.
            name_field: Column used for record names when no schema index is configured.
                Falls back to `name` if absent. If neither exists, records are created
                without names.

        Examples:

            Create a new type and import records::

                records = ln.Record.from_dataframe(df, type="my_df").save()

            Import records into an existing type::

                records = ln.Record.from_dataframe(df, type=sample_sheet).save()

        """
        import pandas as pd

        from .schema import Schema

        if not isinstance(df, pd.DataFrame):
            raise TypeError("`df` needs to be a pandas DataFrame.")
        resolved_type: Record
        if isinstance(type, str):
            imports_type = cls.filter(uid=IMPORTS_UID).one_or_none()
            if imports_type is None:
                imports_type = cls(name="Imports", is_type=True)
                imports_type.uid = IMPORTS_UID
                imports_type = imports_type.save()
            existing_type = cls.filter(
                name=type, is_type=True, type=imports_type
            ).one_or_none()
            if existing_type is not None:
                raise ValueError(
                    f"type '{type}' already exists under 'Imports', please pass it as a Record object to reuse."
                )
            imports_schema = Schema.filter(uid=SCHEMA_IMPORTS_UID).one_or_none()
            if imports_schema is None:
                imports_schema = Schema(name="Imports", is_type=True)
                imports_schema.uid = SCHEMA_IMPORTS_UID
                imports_schema = imports_schema.save()
            inferred_schema = Schema.from_dataframe(df, name=type)
            if inferred_schema is None:
                raise ValueError(
                    "Could not infer a schema from dataframe columns. "
                    "Ensure dataframe columns map to existing Features, or pass an existing Record type object."
                )
            inferred_schema.type = imports_schema
            inferred_schema = inferred_schema.save()
            resolved_type = cls(
                name=type,
                is_type=True,
                type=imports_type,
                schema=inferred_schema,
            ).save()
        else:
            resolved_type = type
        if not resolved_type.is_type:
            raise ValueError("`type` needs to be a record type (`is_type=True`).")
        if resolved_type.name is None:
            raise ValueError("`type` needs to have a non-null `name`.")

        return RecordBatch(
            cls=cls,
            df=df,
            resolved_type=resolved_type,
            name_field=name_field,
        )

    @property
    def features(self) -> FeatureManager:
        """Manage the linked feature values.

        For examples, see :class:`~lamindb.Record`.
        """
        from ._feature_manager import FeatureManager

        return FeatureManager(self)

    @property
    def is_sheet(self) -> bool:
        """Check if record is a `sheet`, i.e., `self.is_type and self.schema is not None`."""
        return self.schema is not None and self.is_type

    def query_parents(self) -> QuerySet:
        """Query all parents of a record recursively.

        While `.parents` retrieves the direct parents, this method
        retrieves all ancestors of the current record.
        """
        return _query_relatives([self], "parents")  # type: ignore

    def query_children(self) -> QuerySet:
        """Query all children of a record recursively.

        While `.children` retrieves the direct children, this method
        retrieves all descendants of a parent.
        """
        return _query_relatives([self], "children")  # type: ignore

    def query_records(self) -> QuerySet:
        """Query records of sub types.

        While `.records` retrieves the records with the current type, this method
        also retrieves sub types and the records with sub types of the current type.
        """
        return _query_relatives([self], "records")  # type: ignore

    def _set_export_run(self, is_run_input: bool | Run | None = None) -> None:
        from lamindb.core._context import context
        from lamindb.models import Run, Transform

        if isinstance(is_run_input, Run):
            run = is_run_input
        elif is_run_input in {True, None}:
            # If there is an active run context, link it as initiator of the
            # internal record export run.
            initiated_by_run = context.run
            # Compatibility path for older instances:
            # historically, "__lamindb_record_export__" transforms could have
            # arbitrary UIDs. We now standardize on a fixed UID to make creation
            # idempotent under concurrency and to avoid duplicate internal
            # transforms. This follows the fixed-UID pattern already used by
            # "save_vitessce_config" (integrations/_vitessce.py) and
            # "__lamindb_transfer__/{instance_uid}" (models/sqlrecord.py).
            # After a few lamindb release cycles (once legacy UIDs are no
            # longer expected), this normalization branch can be removed.
            export_transform_uid = "v6KpQx9mRt2B0000"
            transform = Transform.objects.filter(uid=export_transform_uid).first()
            if transform is None:
                transform = (
                    Transform.objects.filter(
                        key="__lamindb_record_export__", kind="function"
                    )
                    .order_by("created_at")
                    .first()
                )
            if transform is None:
                transform, _ = Transform.objects.get_or_create(
                    uid=export_transform_uid,
                    defaults={
                        "key": "__lamindb_record_export__",
                        "kind": "function",
                    },
                )
            elif transform.uid != export_transform_uid:
                transform.uid = export_transform_uid
                transform.save()
            # Export is treated as a discrete user action, so always create a new
            # run. Transfer reuses runs to avoid repeated sync bookkeeping runs.
            run = Run(
                transform=transform,
                initiated_by_run=initiated_by_run,
                status="started",
            ).save()  # type: ignore
            run.initiated_by_run = initiated_by_run  # available in memory
        else:
            run = None
        self._export_run = run

    @class_and_instance_method
    def to_dataframe(
        cls_or_self,
        recurse: bool = False,
        filters: Any | None = None,
        is_run_input: bool | Run | None = None,
        link_individual_inputs: bool = True,
        **kwargs,
    ) -> pd.DataFrame:
        """Export to a pandas DataFrame.

        This is roughly equivalent to::

            ln.Record.filter(type=sample_type).to_dataframe(include="features")

        `to_dataframe()` ensures that the columns are ordered according to the schema of the type and encodes fields like `uid` and `name`.

        When the sheet schema defines :attr:`~lamindb.Schema.index`, the index feature is
        placed on `df.index` (named after the feature) and encoded metadata columns
        (`__lamindb_record_id__`, `__lamindb_record_uid__`, `__lamindb_record_name__`, etc.)
        are omitted. Sheets without `index` keep the previous export behavior.

        It will also track the record as an input to the current run.

        Example:

            Export all records on a sheet::

                sample_sheet.to_dataframe()

            Export only records with high GC content::

                sample_sheet.to_dataframe(filters=gc_content > 0.55)

        Args:
            recurse: Whether to include records of sub-types recursively.
            filters: Filters applied before export. Supports filter kwargs via
                a `dict`, Django `Q` expressions, and feature predicates
                (e.g. `my_feature > 5`), including iterables of expressions.
            is_run_input: Whether to track the record as a run input.
            link_individual_inputs: Whether to link all exported records as
                inputs of the export run. If `False`, only links the record type.
            **kwargs: Keyword arguments passed to :meth:`~lamindb.models.QuerySet.to_dataframe`.
        """
        import pandas as pd

        if isinstance(cls_or_self, type):
            return type(cls_or_self).to_dataframe(cls_or_self, **kwargs)  # type: ignore
        if not cls_or_self.is_type:
            raise TypeError(
                "to_dataframe() can only be called on the class or on record type instance."
            )
        self = cls_or_self
        assert self.is_type, "Only types can be exported as dataframes"  # noqa: S101

        branch_ids = get_default_branch_ids()
        qs = (
            self.query_records()
            if recurse
            else self.records.filter(branch_id__in=branch_ids)
        )
        if isinstance(filters, dict):
            # Keep kwargs behavior aligned with the historic `self.records.filter(...)`
            # semantics where fields like `name` target record fields.
            qs = qs.filter(**filters)
        elif filters is not None:
            # Feature predicates are only supported through Record.filter(...).
            qs = (
                self.query_records()
                if recurse
                else self.__class__.filter(type=self, branch_id__in=branch_ids)
            )
            if isinstance(filters, Iterable) and not isinstance(filters, (str, bytes)):
                qs = qs.filter(*filters)
            else:
                qs = qs.filter(filters)
        logger.important(f"exporting {qs.count()} records of '{self.name}'")
        if "order_by" not in kwargs:
            kwargs["order_by"] = "id"
        features_arg = kwargs.pop("features", "queryset")
        index_feature = self.schema.index if self.schema is not None else None
        include_record_metadata = export_includes_record_metadata(self.schema)
        df = qs.to_dataframe(
            features=features_arg,
            limit=None,
            record_metadata=include_record_metadata,
            **kwargs,
        )
        if not include_record_metadata and self.schema is not None:
            schema_feature_names = set(
                self.schema.members.values_list("name", flat=True)
            )
            pk_name = self._meta.pk.name
            if pk_name in df.columns and pk_name not in schema_feature_names:
                df = df.drop(columns=[pk_name])
        encoded_id = encode_lamindb_fields_as_columns(self.__class__, "id")
        assert isinstance(encoded_id, str)  # noqa: S101
        encoded_uid = encode_lamindb_fields_as_columns(self.__class__, "uid")
        encoded_name = encode_lamindb_fields_as_columns(self.__class__, "name")
        assert isinstance(encoded_name, str)  # noqa: S101
        # encode the django id, uid and name fields
        if include_record_metadata and df.index.name == "id":
            df.index.name = encoded_id
        if (
            include_record_metadata
            and "uid" in df.columns
            and encoded_uid not in df.columns
        ):
            df = df.rename(columns={"uid": encoded_uid})
        if index_feature is not None:
            if "name" in df.columns and index_feature.name != "name":
                df[index_feature.name] = df["name"]
                df = df.drop(columns=["name"])
            if encoded_name in df.columns:
                df = df.drop(columns=[encoded_name])
            df = apply_schema_index_to_export_dataframe(
                df,
                index_feature,
                encoded_id=encoded_id,
                encoded_name=encoded_name,
                include_record_metadata=include_record_metadata,
            )
        elif "name" in df.columns and encoded_name not in df.columns:
            df = df.rename(columns={"name": encoded_name})
        if not include_record_metadata:
            df = drop_record_metadata_columns(df)
        if self.schema is not None:
            all_features = self.schema.members.all()
            index_feature_uid = None if index_feature is None else index_feature.uid
            desired_order = [
                feature.name
                for feature in all_features
                if index_feature_uid is None or feature.uid != index_feature_uid
            ]
            for feature in all_features:
                if index_feature_uid is not None and feature.uid == index_feature_uid:
                    continue
                if feature.name not in df.columns:
                    df[feature.name] = pd.Series(
                        dtype=convert_to_pandas_dtype(feature._dtype_str)
                    )
        else:
            # sort alphabetically for now
            desired_order = df.columns[2:].tolist()
            desired_order.sort()
        df = reorder_subset_columns_in_df(df, desired_order, position=0)  # type: ignore
        self._set_export_run(is_run_input=is_run_input)
        # always link the type record
        self._export_run.input_records.add(self)
        if link_individual_inputs:
            input_record_ids = qs.values_list("id", flat=True)
            self._export_run.input_records.add(*input_record_ids)
        self._export_run.finished_at = datetime.now(timezone.utc)
        self._export_run._status_code = 0  # completed
        self._export_run.save()
        return df.sort_index()  # order by id

    def to_artifact(
        self,
        key: str | None = None,
        suffix: str | None = None,
        filters: Any | None = None,
        is_run_input: bool | Run | None = None,
        link_individual_inputs: bool = True,
        **kwargs,
    ) -> Artifact:
        """Calls `to_dataframe()` to create an artifact.

        The format defaults to `.csv` unless the key specifies another format or suffix is passed.

        The `key` defaults to `sheet_exports/{self.name}{suffix}` unless a `key` is passed.

        When the sheet schema defines :attr:`~lamindb.Schema.index`, the CSV is written
        with `index=True` so the index feature is preserved on export.

        Example:

            Export all records on a sheet to an artifact::

                sample_sheet.to_artifact()

            Export only records with high GC content::

                sample_sheet.to_artifact(filters=gc_content > 0.55)

        Args:
            key: `str | None = None` The artifact key.
            suffix: `str | None = None` The suffix to append to the default key if no key is passed.
            filters: Filters applied before export.
            is_run_input: Whether to track the record as a run input.
            link_individual_inputs: Whether to link all exported records as
                inputs of the export run. If `False`, only links the record type.
            **kwargs: Keyword arguments passed to :meth:`~lamindb.models.Record.to_dataframe`.
        """
        assert self.is_type, "Only types can be exported as artifacts."
        assert key is None or suffix is None, "Only one of key or suffix can be passed."
        if key is None:
            suffix = ".csv" if suffix is None else suffix
            key = f"sheet_exports/{self.name}{suffix}"
        description = f": {self.description}" if self.description is not None else ""
        return Artifact.from_dataframe(
            self.to_dataframe(
                filters=filters,
                is_run_input=is_run_input,
                link_individual_inputs=link_individual_inputs,
                **kwargs,
            ),
            key=key,
            description=f"Export of sheet {self.uid}{description}",
            schema=self.schema,
            csv_kwargs={
                "index": self.schema is not None and self.schema.index is not None
            },
            run=self._export_run,
        ).save()


# for storing JSON values in records
class RecordJson(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_json")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordjson")
    value: Any = JSONField(default=None, db_default=None)

    class Meta:
        app_label = "lamindb"
        # a list is modeled as a list in json, hence no multi-value association for the same feature unlike for
        # categorical/relational values
        unique_together = ("record", "feature")


# for storing record-like values in records
class RecordRecord(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_record")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrecord")
    value: Record = ForeignKey(Record, PROTECT, related_name="links_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for storing ulabel-like values in records
class RecordULabel(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_ulabel")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordulabel")
    value: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one ulabel per feature, because we likely don't want to have Many
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for storing user-like values in records
class RecordUser(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_user")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recorduser")
    value: User = ForeignKey(User, PROTECT, related_name="links_record")

    class Meta:
        # allows linking exactly one record to one user per feature, because we likely don't want to have Many
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for storing run-like values in records
class RecordRun(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_run")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordrun")
    value: Run = ForeignKey(Run, PROTECT, related_name="links_in_record")

    class Meta:
        # allows linking several records to a single run for the same feature because we'll likely need this
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating runs with records
class RunRecord(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_record")
    record: Record = ForeignKey(Record, PROTECT, related_name="links_run")
    feature: Feature = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_runrecord"
    )
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("run", "record", "feature")


# for storing artifact-like values in records
class RecordArtifact(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_artifact")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="links_recordartifact")
    value: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_in_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating artifacts with records
class ArtifactRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_record")
    record: Record = ForeignKey(Record, PROTECT, related_name="links_artifact")
    feature: Feature = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactrecord"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("artifact", "record", "feature")


# for storing collection-like values in records
class RecordCollection(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_collection")
    feature: Feature = ForeignKey(
        Feature, PROTECT, related_name="links_recordcollection"
    )
    value: Collection = ForeignKey(Collection, PROTECT, related_name="links_in_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating collections with records
class CollectionRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_record"
    )
    record: Record = ForeignKey(Record, PROTECT, related_name="links_collection")
    feature: Feature = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_collectionrecord"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("collection", "record", "feature")


# for storing transform-like values in records
class RecordTransform(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    record: Record = ForeignKey(Record, CASCADE, related_name="values_transform")
    feature: Feature = ForeignKey(
        Feature, PROTECT, related_name="links_recordtransform"
    )
    value: Transform = ForeignKey(Transform, PROTECT, related_name="links_in_record")

    class Meta:
        app_label = "lamindb"
        unique_together = ("record", "feature", "value")


# for annotating transforms with records
class TransformRecord(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_record")
    record: Record = ForeignKey(Record, PROTECT, related_name="links_transform")
    feature: Feature = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_transformrecord"
    )
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now()
    )
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("transform", "record", "feature")
