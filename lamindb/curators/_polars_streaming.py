from __future__ import annotations

from threading import Lock
from typing import TYPE_CHECKING, Literal

import polars as pl
from lamindb_setup import logger

from lamindb.errors import ValidationError

from .core import DataFrameCurator

if TYPE_CHECKING:
    from lamindb.models import Schema


class ValidatedLazyFrame:
    """Validate in-memory sink batches, retaining only annotation metadata."""

    def __init__(
        self,
        frame: pl.LazyFrame,
        schema: Schema,
        features=None,
        on_validation_error: Literal["raise", "keep"] = "raise",
    ):
        self._schema = schema
        self._features = features
        self._members = schema.members.to_list()
        self._lock = Lock()
        self._on_validation_error = on_validation_error
        self.error: str | None = None
        empty = pl.DataFrame(schema=frame.collect_schema())
        self.curator = self._curator(empty)
        if (
            on_validation_error == "keep"
            and self.curator._atomic_curator._coercion_requested
        ):
            raise ValueError(
                "on_validation_error='keep' requires a non-coercing schema: "
                "disable Schema.coerce and Feature.coerce, or use 'raise'."
            )
        self._validate(self.curator)
        self.frame = frame.map_batches(
            self._validate_batch,
            schema=(
                empty.schema
                if on_validation_error == "keep"
                else self._validated_batch(self.curator).schema
            ),
            predicate_pushdown=False,
            projection_pushdown=False,
            slice_pushdown=False,
            streamable=True,
        )

    def _curator(self, batch: pl.DataFrame) -> DataFrameCurator:
        return DataFrameCurator(
            batch, self._schema, features=self._features, members=self._members
        )

    def _validate(self, curator: DataFrameCurator) -> bool:
        try:
            curator.validate()
        except ValidationError as error:
            with self._lock:
                if self.error is None:
                    # Bound diagnostics without retaining exception tracebacks or batches.
                    self.error = str(error)[:4096]
                    if self._on_validation_error == "keep":
                        logger.warning(
                            "Validation failed; continuing the sink without a "
                            f"validating schema. First error: {self.error}"
                        )
            if self._on_validation_error == "raise":
                raise
            return False
        return True

    @staticmethod
    def _validated_batch(curator: DataFrameCurator) -> pl.DataFrame:
        validated = curator.dataset
        if isinstance(validated, pl.LazyFrame):
            validated = validated.collect()
        if not isinstance(validated, pl.DataFrame):
            raise TypeError("The Polars curator must return a Polars dataframe.")
        return validated

    def _validate_batch(self, batch: pl.DataFrame) -> pl.DataFrame:
        from lamindb import settings

        if self._on_validation_error == "keep" and self.error is not None:
            return batch
        curator = self._curator(batch)
        if not self._validate(curator):
            return batch
        validated = (
            batch
            if self._on_validation_error == "keep"
            else self._validated_batch(curator)
        )
        with self._lock:
            for key, vector in curator.cat._cat_vectors.items():
                target = self.curator.cat._cat_vectors[key]
                records = {
                    (record.__class__, record.pk): record
                    for record in target.records or []
                }
                for record in vector.records or []:
                    if len(records) > settings.annotation.n_max_records:
                        break
                    records[(record.__class__, record.pk)] = record
                target.records = list(records.values())
        return validated
