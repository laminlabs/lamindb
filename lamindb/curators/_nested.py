"""Recursive dictionary validation through the existing dataframe curator."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from lamindb.core._compat import is_polars_dataframe
from lamindb.errors import InvalidArgument, ValidationError
from lamindb.models import Schema
from lamindb.models.feature import dtype_schema_uid, validate_nested_schema

if TYPE_CHECKING:
    from lamindb.models import Feature


def validate_nested_column(
    dataset, feature: Feature, *, ancestors: tuple[str, ...] = ()
) -> None:
    """Validate a typed dictionary column without replacing its source values."""
    from .core import ComponentCurator, convert_dict_to_dataframe_for_validation

    schema_uid = dtype_schema_uid(feature._dtype_str)
    schema = Schema.objects.using(feature._state.db).get(uid=schema_uid)
    validate_nested_schema(schema)
    if schema.uid in ancestors:
        raise InvalidArgument(f"Cyclic nested dtype at feature '{feature.name}'.")
    members = schema.members.to_list()
    if feature.coerce or schema.coerce or any(member.coerce for member in members):
        raise InvalidArgument(
            f"Nested feature '{feature.name}' requests coercion. "
            "Nested coercion is not supported; normalize the source values first."
        )
    ancestors = (*ancestors, schema.uid)

    def validate_frame(frame):
        if not members:
            columns = (
                frame.collect_schema().names()
                if is_polars_dataframe(frame)
                else list(frame.columns)
            )
            if schema.maximal_set and columns:
                raise ValidationError(f"Unexpected nested fields: {columns}")
            return
        if is_polars_dataframe(frame):
            import polars as pl

            frame_schema = frame.collect_schema()
            null_checks = []
            for member in members:
                if member.nullable or member.name not in frame_schema:
                    continue
                missing = pl.col(member.name).is_null()
                if frame_schema[member.name].is_float():
                    missing = missing | pl.col(member.name).is_nan()
                null_checks.append(missing.any().alias(member.name))
            if null_checks:
                # Fail before Pandera collects per-row null failure masks.
                failures = frame.select(null_checks).collect().row(0, named=True)
                invalid = [name for name, failed in failures.items() if failed]
                if invalid:
                    raise ValidationError(
                        f"Non-nullable nested fields contain null values: {invalid}"
                    )
        curator = ComponentCurator(
            frame,
            schema,
            using=feature._state.db,
            members=members,
            _nested_ancestors=ancestors,
        )
        curator.validate()

    def validate_pandas(values, dtype_str, nullable):
        present = []
        for value in values:
            if pd.api.types.is_scalar(value) and pd.isna(value):
                if not nullable:
                    raise ValidationError("Null value in a non-nullable nested field.")
            else:
                present.append(value)
        if not present:
            return
        if dtype_str.startswith("list["):
            items = []
            for value in present:
                if not isinstance(value, (list, tuple, np.ndarray, pd.Series)) or (
                    isinstance(value, np.ndarray) and value.ndim == 0
                ):
                    raise ValidationError("Expected a list.")
                items.extend(value)
            validate_pandas(items, dtype_str[5:-1], True)
            return
        rows = []
        for value in present:
            if not isinstance(value, dict):
                raise ValidationError("Expected a dictionary.")
            rows.append(
                convert_dict_to_dataframe_for_validation(
                    value, schema, members=members, index=None
                )
            )
        # A Struct has one field layout, including fields absent in some rows.
        columns = dict.fromkeys(name for row in rows for name in row.columns)
        frame = pd.DataFrame(index=range(len(rows)))
        for name in columns:
            values = [row[name].iloc[0] if name in row else None for row in rows]
            frame[name] = (
                pd.Series(pd.array(values))
                if all(pd.api.types.is_scalar(value) for value in values)
                else pd.Series(values, dtype=object)
            )
        for member in members:
            if member.name not in frame:
                continue
            if member._dtype_str.startswith("cat"):
                frame[member.name] = pd.Categorical(frame[member.name])
            elif member._dtype_str == "dict":
                frame[member.name] = frame[member.name].astype(object)
            if any(row.attrs.get(member.name) == "list_of_categories" for row in rows):
                frame.attrs[member.name] = "list_of_categories"
        validate_frame(frame)

    def validate_polars(frame, dtype_str, nullable):
        import polars as pl

        dtype = frame.collect_schema()[feature.name]
        values = pl.col(feature.name)
        if not nullable and frame.select(values.is_null().any()).collect().item():
            raise ValidationError("Null value in a non-nullable nested field.")
        if dtype == pl.Null:
            return
        if dtype_str.startswith("list["):
            if not isinstance(dtype, pl.List):
                raise ValidationError("Expected a Polars List.")
            items = frame.filter(
                values.is_not_null() & (values.list.len() > 0)
            ).explode(feature.name)
            validate_polars(items, dtype_str[5:-1], True)
            return
        if not isinstance(dtype, pl.Struct):
            raise ValidationError("Expected a Polars Struct.")
        if not frame.select(values.is_not_null().any()).collect().item():
            return
        fields = frame.filter(values.is_not_null()).select(
            [
                values.struct.field(field.name).alias(field.name)
                for field in dtype.fields
            ]
        )
        validate_frame(fields)

    try:
        if is_polars_dataframe(dataset):
            frame = dataset.lazy()
            validate_polars(frame, feature._dtype_str, feature.nullable)
        else:
            validate_pandas(dataset[feature.name], feature._dtype_str, feature.nullable)
    except ValidationError as error:
        raise ValidationError(f"Nested field '{feature.name}': {error}") from error
