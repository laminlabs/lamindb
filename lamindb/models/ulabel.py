from __future__ import annotations

from typing import TYPE_CHECKING, overload

import pgtrigger
from django.conf import settings as django_settings
from django.db import models
from django.db.models import CASCADE, PROTECT

from lamindb.base.fields import (
    CharField,
    DateTimeField,
    ForeignKey,
    TextField,
)
from lamindb.errors import FieldValidationError

from ..base.uids import base62_8
from .can_curate import CanCurate
from .feature import Feature
from .has_parents import HasParents, _query_relatives
from .run import Run, TracksRun, TracksUpdates, User, current_user_id
from .sqlrecord import BaseSQLRecord, HasType, IsLink, SQLRecord, _get_record_kwargs
from .transform import Transform

if TYPE_CHECKING:
    from datetime import datetime

    from .artifact import Artifact
    from .block import ULabelBlock
    from .collection import Collection
    from .project import Project
    from .query_manager import RelatedManager
    from .query_set import QuerySet
    from .record import Record


class ULabel(SQLRecord, HasType, HasParents, CanCurate, TracksRun, TracksUpdates):
    """Universal labels.

    In some cases you may just want to create a simple label, then `ULabel` is for you.

    It behaves like `Record`, just without the ability to store features.

    Args:
        name: `str` A name.
        description: `str | None = None` A description.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement (e.g. column of a sheet, attribute of a record).
        :meth:`~lamindb.ULabel`
            Like `ULabel`, but with the ability to store features.

    Examples:

        Create a label::

            train_split = ln.ULabel(name="train").save()

        Organize ulabels in a hierarchy::

            split_type = ln.ULabel(name="Split", is_type=True).save()
            train_split = ln.ULabel(name="train", type="split_type").save()

        Label an artifact::

            artifact.ulabels.add(train_split)

        Query artifacts by label::

            ln.Artifact.filter(ulabels=train_split).to_dataframe()
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
                    name="prevent_ulabel_type_cycle",
                    operation=pgtrigger.Update | pgtrigger.Insert,
                    when=pgtrigger.Before,
                    condition=pgtrigger.Condition("NEW.type_id IS NOT NULL"),
                    func="""
                        -- Check for direct self-reference
                        IF NEW.type_id = NEW.id THEN
                            RAISE EXCEPTION 'Cannot set type: ulabel cannot be its own type';
                        END IF;

                        -- Check for cycles in the type chain
                        IF EXISTS (
                            WITH RECURSIVE type_chain AS (
                                SELECT type_id, 1 as depth
                                FROM lamindb_ulabel
                                WHERE id = NEW.type_id

                                UNION ALL

                                SELECT r.type_id, tc.depth + 1
                                FROM lamindb_ulabel r
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
        editable=False, unique=True, db_index=True, max_length=8, default=base62_8
    )
    """A universal random id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True)
    """Name or title of ulabel."""
    type: ULabel | None = ForeignKey("self", PROTECT, null=True, related_name="ulabels")
    """Type of ulabel, e.g., `"donor"`, `"split"`, etc. ← :attr:`~lamindb.ULabel.ulabels`

    Allows to group ulabels by type, e.g., all donors, all split ulabels, etc.
    """
    ulabels: ULabel
    """ULabels of this type (can only be non-empty if `is_type` is `True`)."""
    description: str | None = TextField(null=True)
    """A description."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A simple reference like URL or external ID."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of simple reference."""
    parents: RelatedManager[ULabel] = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent entities of this ulabel ← :attr:`~lamindb.ULabel.children`.

    For advanced use cases, you can build an ontology under a given `type`.

    Say, if you modeled `CellType` as a `ULabel`, you would introduce a type `CellType` and model the hiearchy of cell types under it.
    """
    children: RelatedManager[ULabel]
    """Child entities of this ulabel.

    Reverse accessor for parents.
    """
    transforms: RelatedManager[Transform]
    """The transforms annotated by this ulabel ← :attr:`~lamindb.Transform.ulabels`."""
    runs: RelatedManager[Run]
    """The runs annotated by this ulabel ← :attr:`~lamindb.Run.ulabels`."""
    artifacts: RelatedManager[Artifact] = models.ManyToManyField(
        "Artifact", through="ArtifactULabel", related_name="ulabels"
    )
    """The artifacts annotated by this ulabel ← :attr:`~lamindb.Artifact.ulabels`."""
    collections: RelatedManager[Collection]
    """The collections annotated by this ulabel ← :attr:`~lamindb.Collection.ulabels`."""
    projects: RelatedManager[Project]
    """The projects annotating this ulabel ← :attr:`~lamindb.Project.ulabels`."""
    linked_in_records: RelatedManager[Record] = models.ManyToManyField(
        "Record",
        through="RecordULabel",
        related_name="linked_ulabels",
    )
    """Records linking this ulabel as a value ← :attr:`~lamindb.Record.linked_ulabels`."""
    ablocks: ULabelBlock
    """Attached blocks ← :attr:`~lamindb.ULabelBlock.ulabel`."""

    @overload
    def __init__(
        self,
        name: str,
        type: ULabel | None = None,
        is_type: bool = False,
        description: str | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
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
        description: str | None = kwargs.pop("description", None)
        reference: str | None = kwargs.pop("reference", None)
        reference_type: str | None = kwargs.pop("reference_type", None)
        branch = kwargs.pop("branch", None)
        branch_id = kwargs.pop("branch_id", 1)
        space = kwargs.pop("space", None)
        space_id = kwargs.pop("space_id", 1)
        _skip_validation = kwargs.pop("_skip_validation", False)
        _aux = kwargs.pop("_aux", None)
        if len(kwargs) > 0:
            valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(ULabel)])
            raise FieldValidationError(
                f"Only {valid_keywords} are valid keyword arguments"
            )
        super().__init__(
            name=name,
            type=type,
            is_type=is_type,
            description=description,
            reference=reference,
            reference_type=reference_type,
            branch=branch,
            branch_id=branch_id,
            space=space,
            space_id=space_id,
            _skip_validation=_skip_validation,
            _aux=_aux,
        )

    def query_ulabels(self) -> QuerySet:
        """Query ulabels of sub types.

        While `.ulabels` retrieves the ulabels with the current type, this method
        also retrieves sub types and the ulabels with sub types of the current type.
        """
        return _query_relatives([self], "ulabels")  # type: ignore


class ArtifactULabel(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey("Artifact", CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactulabel", default=None
    )

    class Meta:
        # can have the same label linked to the same artifact if the feature is
        # different
        app_label = "lamindb"
        unique_together = ("artifact", "ulabel", "feature")


class TransformULabel(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_transform")

    class Meta:
        app_label = "lamindb"
        unique_together = ("transform", "ulabel")


class RunULabel(BaseSQLRecord, IsLink):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_run")
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""

    class Meta:
        app_label = "lamindb"
        unique_together = ("run", "ulabel")


class CollectionULabel(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        "Collection", CASCADE, related_name="links_ulabel"
    )
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_collection")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_collectionulabel", default=None
    )

    class Meta:
        app_label = "lamindb"
        unique_together = ("collection", "ulabel")
