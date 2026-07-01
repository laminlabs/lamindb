from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Literal, cast

from django.db.models import Case, Q, TextField, Value, When
from django.db.models.functions import Concat
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from upath import UPath

from .artifact import Artifact, track_run_input
from .collection import Collection, _load_concat_artifacts

if TYPE_CHECKING:
    from anndata import AnnData
    from lamindb_setup.types import AnyPathStr
    from pandas import DataFrame
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset

    from ..core._mapped_collection import MappedCollection
    from .record import Record
    from .run import Run


UNORDERED_WARNING = (
    "this query set is unordered, consider using `.order_by()` first "
    "to avoid opening the artifacts in an arbitrary order"
)


# maybe make this abstract
class ArtifactSet(Iterable):
    """Abstract class representing sets of artifacts returned by queries.

    This class automatically extends :class:`~lamindb.models.BasicQuerySet`
    and :class:`~lamindb.models.QuerySet` when the base model is :class:`~lamindb.Artifact`.

    Examples:

        >>> artifacts = ln.Artifact.filter(otype="AnnData")
        >>> artifacts # an instance of ArtifactQuerySet inheriting from ArtifactSet
    """

    @doc_args(Collection.load.__doc__)
    def load(
        self,
        join: Literal["inner", "outer"] = "outer",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> DataFrame | AnnData:
        """{}"""  # noqa: D415
        if not self.ordered:  # type: ignore
            logger.warning(UNORDERED_WARNING)

        artifacts: list[Artifact] = list(self)
        concat_object = _load_concat_artifacts(artifacts, join, **kwargs)
        # track only if successful
        track_run_input(artifacts, is_run_input)
        return concat_object

    @doc_args(Collection.open.__doc__)
    def open(
        self,
        engine: Literal["pyarrow", "polars"] = "pyarrow",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> PyArrowDataset | Iterator[PolarsLazyFrame]:
        """{}"""  # noqa: D415
        from ..core.storage._backed_access import _open_dataframe

        if not self.ordered:  # type: ignore
            logger.warning(UNORDERED_WARNING)

        artifacts: list[Artifact] = list(self)
        paths: list[UPath] = [artifact.path for artifact in artifacts]

        dataframe = _open_dataframe(paths, engine=engine, **kwargs)
        # track only if successful
        track_run_input(artifacts, is_run_input)
        return dataframe

    @doc_args(Collection.mapped.__doc__)
    def mapped(
        self,
        layers_keys: str | list[str] | None = None,
        obs_keys: str | list[str] | None = None,
        obsm_keys: str | list[str] | None = None,
        obs_filter: dict[str, str | list[str]] | None = None,
        join: Literal["inner", "outer"] | None = "inner",
        encode_labels: bool | list[str] = True,
        unknown_label: str | dict[str, str] | None = None,
        cache_categories: bool = True,
        parallel: bool = False,
        dtype: str | None = None,
        stream: bool = False,
        is_run_input: bool | None = None,
    ) -> MappedCollection:
        """{}"""  # noqa: D415
        from ..core._mapped_collection import MappedCollection

        if not self.ordered:  # type: ignore
            logger.warning(UNORDERED_WARNING)

        artifacts: list[Artifact] = []
        paths: list[UPath] = []
        for artifact in self:
            if ".h5ad" not in artifact.suffix and ".zarr" not in artifact.suffix:
                logger.warning(f"ignoring artifact with suffix {artifact.suffix}")
                continue
            elif not stream:
                paths.append(artifact.cache())
            else:
                paths.append(artifact.path)
            artifacts.append(artifact)
        ds = MappedCollection(
            paths,
            layers_keys,
            obs_keys,
            obsm_keys,
            obs_filter,
            join,
            encode_labels,
            unknown_label,
            cache_categories,
            parallel,
            dtype,
        )
        # track only if successful
        track_run_input(artifacts, is_run_input)
        return ds


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
        features: bool | str | list[str] | None = None,
        limit: int | None = 20,
        order_by: str | None = "-id",
        record_metadata: bool = True,
        is_run_input: bool | Run | None = None,
        link_individual_inputs: bool = True,
        _record_type: Record | None = None,
    ) -> DataFrame:
        import pandas as pd

        from .feature import convert_to_pandas_dtype
        from .query_set import (
            SEARCH_QUERY_DEFAULT_LIMIT,
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
                return BasicQuerySet.to_dataframe(
                    qs,
                    include=include,
                    features=features,
                    limit=limit,
                    order_by=order_by,
                    record_metadata=record_metadata,
                )

            # `type_id` points to record types (`is_type=True`) by model design.
            record_type = Record.get(id=type_ids[0])

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

        record_type._set_export_run(is_run_input=is_run_input)
        export_run = record_type._export_run
        if export_run is not None:
            export_run.input_records.add(record_type)
            if link_individual_inputs:
                input_record_ids = qs.values_list("id", flat=True)
                export_run.input_records.add(*input_record_ids)
            from datetime import datetime, timezone

            export_run.finished_at = datetime.now(timezone.utc)
            export_run._status_code = 0
            export_run.save()
        qs._record_export_run = export_run
        return df.sort_index()


def artifacts_from_path(artifacts: ArtifactSet, path: AnyPathStr) -> ArtifactSet:
    """Returns artifacts in the query set that are registered for the provided path."""
    from lamindb.models import BasicQuerySet, QuerySet

    # not QuerySet but only BasicQuerySet
    assert isinstance(artifacts, BasicQuerySet) and not isinstance(artifacts, QuerySet)  # noqa: S101

    upath = UPath(path)

    path_str = upath.as_posix()

    stem = upath.stem
    stem_len = len(stem)
    suffix = upath.suffix

    if stem_len == 16:
        qs = artifacts.filter(
            Q(_key_is_virtual=True) | Q(key__isnull=True),
            _real_key__isnull=True,
            _overwrite_versions=True,
            suffix=suffix,
            uid__startswith=stem,
        )
    elif stem_len == 20:
        qs = artifacts.filter(
            Q(_key_is_virtual=True) | Q(key__isnull=True),
            _real_key__isnull=True,
            _overwrite_versions=False,
            suffix=suffix,
            uid=stem,
        )
    else:
        qs = None

    if qs:  # an empty query set evaluates to False
        return qs

    qs = (
        artifacts.filter(Q(_key_is_virtual=False) | Q(_real_key__isnull=False))
        .alias(
            db_path=Case(
                When(
                    _real_key__isnull=False,
                    then=Concat(
                        "storage__root",
                        Value("/"),
                        "_real_key",
                        output_field=TextField(),
                    ),
                ),
                default=Concat(
                    "storage__root", Value("/"), "key", output_field=TextField()
                ),
                output_field=TextField(),
            )
        )
        .filter(db_path=path_str)
    )

    return qs
