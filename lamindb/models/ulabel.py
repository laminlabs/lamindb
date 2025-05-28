from __future__ import annotations

from typing import TYPE_CHECKING, overload

from django.db import models
from django.db.models import CASCADE, PROTECT

from lamindb.base import deprecated
from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)
from lamindb.errors import FieldValidationError

from ..base.ids import base62_8
from .can_curate import CanCurate
from .feature import Feature
from .has_parents import HasParents
from .run import Run, TracksRun, TracksUpdates, User, current_user_id
from .sqlrecord import BaseSQLRecord, IsLink, SQLRecord, _get_record_kwargs
from .transform import Transform

if TYPE_CHECKING:
    from datetime import datetime

    from .artifact import Artifact
    from .collection import Collection
    from .project import Project


class ULabel(SQLRecord, HasParents, CanCurate, TracksRun, TracksUpdates):
    """Universal labels.

    Args:
        name: `str` A name.
        description: `str` A description.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    A `ULabel` record provides the easiest way to annotate a dataset
    with a label: `"My project"`, `"curated"`, or `"Batch X"`:

        >>> my_project = ULabel(name="My project").save()
        >>> artifact.ulabels.add(my_project)

    Often, a ulabel is measured *within* a dataset. For instance, an artifact
    might characterize 2 species of the Iris flower (`"setosa"` &
    `"versicolor"`) measured by a `"species"` feature. Use the
    :class:`~lamindb.curators.DataFrameCurator` flow to automatically parse, validate, and
    annotate with labels that are contained in `DataFrame` objects.

    .. note::

        If you work with complex entities like cell lines, cell types, tissues,
        etc., consider using the pre-defined biological registries in
        :mod:`bionty` to label artifacts & collections.

        If you work with biological samples, likely, the only sustainable way of
        tracking metadata, is to create a custom schema module.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement for artifacts & collections.
        :attr:`~lamindb.Artifact.features`
            Feature manager for an artifact.

    Examples:

        Create a ulabel:

        >>> train_split = ln.ULabel(name="train").save()

        Organize ulabels in a hierarchy:

        >>> split_type = ln.ULabel(name="Split", is_type=True).save()
        >>> train_split = ln.ULabel(name="train", type="split_type").save()

        Label an artifact:

        >>> artifact.ulabels.add(ulabel)

        Query an artifact by ulabel:

        >>> ln.Artifact.filter(ulabels=train_split).df()
    """

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

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
    """Type of ulabel, e.g., `"donor"`, `"split"`, etc.

    Allows to group ulabels by type, e.g., all donors, all split ulabels, etc.
    """
    ulabels: ULabel
    """ULabels of this type (can only be non-empty if `is_type` is `True`)."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type.

    For example, a ulabel "Project" would be a type, and the actual projects "Project 1", "Project 2", would be records of that `type`.
    """
    description: str | None = CharField(null=True, db_index=True)
    """A description (optional)."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A simple reference like URL or external ID."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of simple reference."""
    parents: ULabel = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent entities of this ulabel.

    For advanced use cases, you can build an ontology under a given `type`.

    Say, if you modeled `CellType` as a `ULabel`, you would introduce a type `CellType` and model the hiearchy of cell types under it.
    """
    children: ULabel
    """Child entities of this ulabel.

    Reverse accessor for parents.
    """
    transforms: Transform
    """Linked transforms."""
    runs: Run
    """Linked runs."""
    artifacts: Artifact
    """Linked artifacts."""
    collections: Collection
    """Linked collections."""
    projects: Project
    """Linked projects."""

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
            _skip_validation=_skip_validation,
            _aux=_aux,
        )

    @property
    @deprecated("ulabels")
    def records(self) -> list[ULabel]:
        """Return all instances of this type."""
        return self.ulabels


class ArtifactULabel(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey("Artifact", CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactulabel", default=None
    )
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        # can have the same label linked to the same artifact if the feature is
        # different
        unique_together = ("artifact", "ulabel", "feature")


class TransformULabel(BaseSQLRecord, IsLink, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_transform")

    class Meta:
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
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        unique_together = ("collection", "ulabel")
