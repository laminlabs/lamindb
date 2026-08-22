from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, cast

from lamin_utils import logger

from .artifact import Artifact
from .query_manager import SEARCH_QUERY_DEFAULT_LIMIT

if TYPE_CHECKING:
    from pandas import DataFrame

    from .record import Record
    from .run import Run


class RecordSet(Iterable):
    """Abstract class representing sets of records returned by queries.

    This class automatically extends :class:`~lamindb.models.BasicQuerySet`
    and :class:`~lamindb.models.QuerySet` when the base model is
    :class:`~lamindb.Record`.
    """

    def to_dataframe(
        self,
        *,
        include: str | list[str] | None = None,
        features: str | list[str] | None = None,
        limit: int | None = SEARCH_QUERY_DEFAULT_LIMIT,
        order_by: str | None = "-id",
        record_metadata: bool = True,
        is_run_input: bool | Run | None = None,
        link_individual_inputs: bool = True,
        use_export_run: bool = False,
        _record_type: Record | None = None,
    ) -> DataFrame:
        """Export records in the queryset to a pandas DataFrame.

        Run-input linking is only performed when feature export is requested,
        i.e., `include="features"` (or `include` contains `"features"`).
        If features are not included, the method falls back to the generic
        queryset export path and does not link inputs.

        Args:
            include: Fields to include. Pass `"features"` (or include it in a
              list) for sheet-style export and optional run-input linking.
            features: Feature names to include when exporting features.
            limit: Maximum number of records to export.
            order_by: Ordering for exported records.
            record_metadata: Whether to include encoded record metadata columns.
            is_run_input: Whether to link exported records as run inputs.
            link_individual_inputs: Whether to link each exported record as an input.
            use_export_run: Whether to use a dedicated internal export run.
        """
        import pandas as pd

        from .feature import convert_to_pandas_dtype
        from .query_set import (
            BasicQuerySet,
            encode_lamindb_fields_as_columns,
            reorder_subset_columns_in_df,
        )
        from .record import (
            Record,
            apply_schema_index_to_export_dataframe,
            drop_record_metadata_columns,
            export_includes_record_metadata,
        )

        qs = cast(BasicQuerySet, self)
        qs._record_export_type = None
        qs._record_export_run = None

        include_features = include == "features" or (
            isinstance(include, list) and "features" in include
        )
        if not include_features:
            return BasicQuerySet.to_dataframe(
                qs,
                include=include,
                features=features,
                limit=limit,
                order_by=order_by,
                record_metadata=record_metadata,
            )

        if _record_type is not None:
            record_type = _record_type
        else:
            type_ids = list(qs.values_list("type_id", flat=True).distinct()[:2])
            if len(type_ids) != 1 or type_ids[0] is None:
                # Sheet-style export requires a single concrete record type context.
                logger.warning(
                    "falling back to generic Record queryset export because the queryset "
                    "does not resolve to exactly one non-null `type_id`"
                )
                qs._record_export_type = None
                qs._record_export_run = None
                return BasicQuerySet.to_dataframe(
                    qs,
                    include=include,
                    features=features,
                    limit=limit,
                    order_by=order_by,
                    record_metadata=record_metadata,
                )

            # `type_id` points to record types (`is_type=True`) by model design.
            # `type_ids` were read from `qs.db`; resolve the type on the same
            # instance so cross-instance reads don't hit the default connection.
            record_type = (
                Record.get(id=type_ids[0])
                if qs.db in (None, "default")
                else Record.connect(qs.db).get(id=type_ids[0])
            )
        qs._record_export_type = record_type

        logger.important(f"exporting {qs.count()} records of '{record_type.name}'")

        features_arg = "queryset" if features is None else features
        limit_arg = None if limit == SEARCH_QUERY_DEFAULT_LIMIT else limit
        order_by_arg = "id" if order_by == "-id" else order_by
        index_feature = (
            record_type.schema.index if record_type.schema is not None else None
        )
        include_record_metadata = export_includes_record_metadata(record_type.schema)
        df = BasicQuerySet.to_dataframe(
            qs,
            include=include,
            features=features_arg,
            limit=limit_arg,
            order_by=order_by_arg,
            record_metadata=include_record_metadata,
        )
        if not include_record_metadata and record_type.schema is not None:
            schema_feature_names = set(
                record_type.schema.members.values_list("name", flat=True)
            )
            pk_name = record_type._meta.pk.name
            if pk_name in df.columns and pk_name not in schema_feature_names:
                df = df.drop(columns=[pk_name])
        encoded_id = encode_lamindb_fields_as_columns(record_type.__class__, "id")
        assert isinstance(encoded_id, str)  # noqa: S101
        encoded_uid = encode_lamindb_fields_as_columns(record_type.__class__, "uid")
        encoded_name = encode_lamindb_fields_as_columns(record_type.__class__, "name")
        assert isinstance(encoded_name, str)  # noqa: S101
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
        if record_type.schema is not None:
            all_features = record_type.schema.members.all()
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
            desired_order = df.columns[2:].tolist()
            desired_order.sort()
        df = reorder_subset_columns_in_df(df, desired_order, position=0)  # type: ignore

        record_type._set_export_run(
            is_run_input=is_run_input,
            use_export_run=use_export_run,
        )
        run_for_input_linking = record_type._export_run
        if run_for_input_linking is not None:
            # avoid cycles: a record cannot be both input and output of the same run
            if record_type.run_id != run_for_input_linking.id:
                run_for_input_linking.input_records.add(record_type)
            if link_individual_inputs:
                input_record_ids = list(
                    qs.exclude(run_id=run_for_input_linking.id).values_list(
                        "id", flat=True
                    )
                )
                if input_record_ids:
                    run_for_input_linking.input_records.add(*input_record_ids)
        if use_export_run and run_for_input_linking is not None:
            from datetime import datetime, timezone

            run_for_input_linking.finished_at = datetime.now(timezone.utc)
            run_for_input_linking._status_code = 0
            run_for_input_linking.save()
        qs._record_export_run = run_for_input_linking if use_export_run else None
        return df.sort_index()

    def to_artifact(
        self,
        key: str | None = None,
        suffix: str | None = None,
        is_run_input: bool | Run | None = None,
        link_individual_inputs: bool = True,
        **kwargs,
    ) -> Artifact:
        """Calls `to_dataframe()` to create an artifact.

        The format defaults to `.csv` unless `suffix` is passed or `key` specifies another format.

        The `key` defaults to `record_exports/{type_name}_subset{suffix}` unless a
        `key` is passed.

        Args:
            key: The artifact key.
            suffix: The suffix to append to the default key if no key is passed.
            is_run_input: Whether to track records as run inputs.
            link_individual_inputs: Whether to link all exported records as
                inputs of the export run. If `False`, only links the record type.
            **kwargs: Keyword arguments passed to :meth:`~lamindb.models.RecordSet.to_dataframe`.
        """
        from .query_set import BasicQuerySet

        assert key is None or suffix is None, "Only one of key or suffix can be passed."
        qs = cast(BasicQuerySet, self)
        kwargs.setdefault("include", "features")
        df = self.to_dataframe(
            is_run_input=is_run_input,
            link_individual_inputs=link_individual_inputs,
            use_export_run=True,
            **kwargs,
        )
        record_type = getattr(qs, "_record_export_type", None)
        export_run = getattr(qs, "_record_export_run", None)
        if export_run is None and record_type is not None:
            export_run = getattr(record_type, "_export_run", None)
        if key is None:
            suffix = ".csv" if suffix is None else suffix
            type_name = record_type.name if record_type is not None else "record"
            key = f"record_exports/{type_name}_subset{suffix}"
        schema = record_type.schema if record_type is not None else None
        return Artifact.from_dataframe(
            df,
            key=key,
            description=(
                f"Export of {record_type.name} subset"
                if record_type is not None
                else "Export of record subset"
            ),
            schema=schema,
            csv_kwargs={"index": schema is not None and schema.index is not None},
            run=export_run,
            space=record_type.space if record_type is not None else None,
        ).save()
