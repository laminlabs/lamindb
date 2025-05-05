from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Iterator

    from anndata import AnnData
    from pandas import DataFrame
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset


class ArtifactQuerySet:
    def load(
        self,
        join: Literal["inner", "outer"] = "outer",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> DataFrame | AnnData:
        pass

    def open(
        self,
        engine: Literal["pyarrow", "polars"] = "pyarrow",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> PyArrowDataset | Iterator[PolarsLazyFrame]:
        pass
