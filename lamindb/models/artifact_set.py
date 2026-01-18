from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Literal

from django.db.models import Case, Q, TextField, Value, When
from django.db.models.functions import Concat
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from upath import UPath

from ..core._mapped_collection import MappedCollection
from ..core.storage._backed_access import _open_dataframe
from .artifact import Artifact, track_run_input
from .collection import Collection, _load_concat_artifacts

if TYPE_CHECKING:
    from anndata import AnnData
    from lamindb_setup.types import UPathStr
    from pandas import DataFrame
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset


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


def artifacts_from_path(artifacts: ArtifactSet, path: UPathStr) -> ArtifactSet:
    """Returns artifacts in the query set that are registered for the provided path."""
    from lamindb.models import BasicQuerySet, QuerySet

    # not QuerySet but only BasicQuerySet
    assert isinstance(artifacts, BasicQuerySet) and not isinstance(artifacts, QuerySet)  # noqa: S101

    upath = UPath(path)

    path_str = upath.as_posix()

    stem = upath.stem
    stem_len = len(stem)

    if stem_len == 16:
        qs = artifacts.filter(
            Q(_key_is_virtual=True) | Q(key__isnull=True),
            _real_key__isnull=True,
            uid__startswith=stem,
        )
    elif stem_len == 20:
        qs = artifacts.filter(
            Q(_key_is_virtual=True) | Q(key__isnull=True),
            _real_key__isnull=True,
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
