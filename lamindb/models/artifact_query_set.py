from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Literal

from lamin_utils import logger

from ..core._mapped_collection import MappedCollection
from ..core.storage._backed_access import _open_dataframe

if TYPE_CHECKING:
    from anndata import AnnData
    from django.models import QuerySet
    from pandas import DataFrame
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset
    from upath import UPath


def _check_ordered_artifacts(qs: QuerySet):
    from ._artifact import Artifact

    if qs.model != Artifact:
        raise ValueError("A query set should consist of artifacts to be opened.")
    if not qs.ordered:
        logger.warning(
            "this query set is unordered, consider using `.order_by()` first "
            "to avoid opening the artifacts in an arbitrary order"
        )


class ArtifactQuerySet(Iterable):
    def load(
        self,
        join: Literal["inner", "outer"] = "outer",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> DataFrame | AnnData:
        from .artifact import Artifact, _track_run_input
        from .collection import _load_concat_artifacts

        _check_ordered_artifacts(self)

        artifacts: list[Artifact] = list(self)
        concat_object = _load_concat_artifacts(artifacts, join, **kwargs)
        # track only if successful
        _track_run_input(artifacts, is_run_input)
        return concat_object

    def open(
        self,
        engine: Literal["pyarrow", "polars"] = "pyarrow",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> PyArrowDataset | Iterator[PolarsLazyFrame]:
        from ._artifact import Artifact, _track_run_input

        _check_ordered_artifacts(self)

        artifacts: list[Artifact] = list(self)
        paths: list[UPath] = [artifact.path for artifact in artifacts]

        dataframe = _open_dataframe(paths, engine=engine, **kwargs)
        # track only if successful
        _track_run_input(artifacts, is_run_input)
        return dataframe

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
        from ._artifact import Artifact, _track_run_input

        _check_ordered_artifacts(self)

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
        _track_run_input(artifacts, is_run_input)
        return ds
