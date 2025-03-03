from datetime import datetime
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Optional,
    overload,
)

from django.db import models
from django.db.models import (
    CASCADE,
    PROTECT,
    Q,
)
from django.db.utils import IntegrityError
from lamindb_setup.core.hashing import HASH_LENGTH, hash_dict
from upath import UPath

from lamindb.base.fields import (
    BooleanField,
    CharField,
    DateTimeField,
    ForeignKey,
)

from ..base.ids import base62_12
from ..errors import (
    ValidationError,
)
from .base import TracksRun, TracksUpdates, current_user_id
from .can_curate import CanCurate
from .record import BasicRecord, Record

if TYPE_CHECKING:
    from .artifact import Artifact
    from .run import Run
    from .schema import Schema
    from .transform import Transform


class User(BasicRecord, CanCurate):
    """Users.

    All data in this registry is synched from `lamin.ai` to ensure a universal
    user identity. There is no need to manually create records.

    Examples:

        Query a user by handle:

        >>> user = ln.User.get(handle="testuser1")
        >>> user
    """

    _name_field: str = "handle"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(editable=False, unique=True, db_index=True, max_length=8)
    """Universal id, valid across DB instances."""
    handle: str = CharField(max_length=30, unique=True, db_index=True)
    """Universal handle, valid across DB instances (required)."""
    name: str | None = CharField(max_length=150, db_index=True, null=True)
    """Name (optional)."""  # has to match hub specification, where it's also optional
    created_artifacts: "Artifact"
    """Artifacts created by user."""
    created_transforms: "Transform"
    """Transforms created by user."""
    created_runs: "Run"
    """Runs created by user."""
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    updated_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of last update to record."""

    @overload
    def __init__(
        self,
        handle: str,
        email: str,
        name: str | None,
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
        super().__init__(*args, **kwargs)


class Storage(Record, TracksRun, TracksUpdates):
    """Storage locations.

    A storage location is either a directory/folder (local or in the cloud) or
    an entire S3/GCP bucket.

    A LaminDB instance can manage and link multiple storage locations. But any
    storage location is managed by *at most one* LaminDB instance.

    .. dropdown:: Managed vs. linked storage locations

        The LaminDB instance can update & delete artifacts in managed storage
        locations but merely read artifacts in linked storage locations.

        When you transfer artifacts from another instance, the default is to
        only copy metadata into the target instance, but merely link the data.

        The `instance_uid` field indicates the managing LaminDB instance of a
        storage location.

        When you delete a LaminDB instance, you'll be warned about data in managed
        storage locations while data in linked storage locations is ignored.

    See Also:
        :attr:`~lamindb.core.Settings.storage`
            Default storage.
        :attr:`~lamindb.setup.core.StorageSettings`
            Storage settings.

    Examples:

        Configure the default storage location upon initiation of a LaminDB instance::

            lamin init --storage ./mydata # or "s3://my-bucket" or "gs://my-bucket"

        View the default storage location:

        >>> ln.settings.storage
        PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata')

        Dynamically change the default storage:

        >>> ln.settings.storage = "./storage_2" # or a cloud bucket
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "root"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, default=base62_12, db_index=True
    )
    """Universal id, valid across DB instances."""
    # we are very conservative here with 255 characters
    root: str = CharField(db_index=True, unique=True)
    """Root path of storage. n s3 path.  local path, etc. (required)."""
    description: str | None = CharField(db_index=True, null=True)
    """A description of what the storage location is used for (optional)."""
    type: str = CharField(max_length=30, db_index=True)
    """Can be "local" vs. "s3" vs. "gs"."""
    region: str | None = CharField(max_length=64, db_index=True, null=True)
    """Cloud storage region, if applicable."""
    instance_uid: str | None = CharField(max_length=12, db_index=True, null=True)
    """Instance that manages this storage location."""
    artifacts: "Artifact"
    """Artifacts contained in this storage location."""

    @overload
    def __init__(
        self,
        root: str,
        type: str,
        region: str | None,
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
        super().__init__(*args, **kwargs)

    @property
    def path(self) -> Path | UPath:
        """Bucket or folder path.

        Cloud storage bucket:

        >>> ln.Storage("s3://my-bucket").save()

        Directory/folder in cloud storage:

        >>> ln.Storage("s3://my-bucket/my-directory").save()

        Local directory/folder:

        >>> ln.Storage("./my-directory").save()
        """
        from lamindb_setup.core.upath import create_path

        access_token = self._access_token if hasattr(self, "_access_token") else None
        return create_path(self.root, access_token=access_token)


class Param(Record, CanCurate, TracksRun, TracksUpdates):
    """Parameters of runs & models."""

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    name: str = CharField(max_length=100, db_index=True)
    dtype: str | None = CharField(db_index=True, null=True)
    """Data type ("num", "cat", "int", "float", "bool", "datetime").

    For categorical types, can define from which registry values are
    sampled, e.g., `cat[ULabel]` or `cat[bionty.CellType]`.
    """
    type: Optional["Param"] = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of param (e.g., 'Pipeline', 'ModelTraining', 'PostProcessing').

    Allows to group features by type, e.g., all read outs, all metrics, etc.
    """
    records: "Param"
    """Records of this type."""
    is_type: bool = BooleanField(default=False, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    _expect_many: bool = models.BooleanField(default=False, db_default=False)
    """Indicates whether values for this param are expected to occur a single or multiple times for an artifact/run (default `False`).

    - if it's `False` (default), the values mean artifact/run-level values and a dtype of `datetime` means `datetime`
    - if it's `True`, the values are from an aggregation, which this seems like an edge case but when characterizing a model ensemble trained with different parameters it could be relevant
    """
    schemas: "Schema" = models.ManyToManyField(
        "Schema", through="SchemaParam", related_name="params"
    )
    """Feature sets linked to this feature."""
    # backward fields
    values: "ParamValue"
    """Values for this parameter."""

    def __init__(self, *args, **kwargs):
        from .feature import process_init_feature_param

        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args, **kwargs)
            return None

        dtype = kwargs.get("dtype", None)
        kwargs = process_init_feature_param(args, kwargs, is_param=True)
        super().__init__(*args, **kwargs)
        dtype_str = kwargs.pop("dtype", None)
        if not self._state.adding:
            if not (
                self.dtype.startswith("cat")
                if dtype == "cat"
                else self.dtype == dtype_str
            ):
                raise ValidationError(
                    f"Feature {self.name} already exists with dtype {self.dtype}, you passed {dtype_str}"
                )


# FeatureValue behaves in many ways like a link in a LinkORM
# in particular, we don't want a _public field on it
# Also, we don't inherit from TracksRun because a ParamValue
# is typically created before a run is created and we want to
# avoid delete cycles (for Model params though it might be helpful)
class ParamValue(Record):
    """Parameter values.

    Is largely analogous to `FeatureValue`.
    """

    # we do not have a unique constraint on param & value because it leads to hashing errors
    # for large dictionaries: https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0000
    # we do not hash values because we have `get_or_create` logic all over the place
    # and also for checking whether the (param, value) combination exists
    # there does not seem an issue with querying for a dict-like value
    # https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0001
    _name_field: str = "value"

    param: Param = ForeignKey(Param, CASCADE, related_name="values")
    """The dimension metadata."""
    value: Any = (
        models.JSONField()
    )  # stores float, integer, boolean, datetime or dictionaries
    """The JSON-like value."""
    # it'd be confusing and hard to populate a run here because these
    # values are typically created upon creating a run
    # hence, ParamValue does _not_ inherit from TracksRun but manually
    # adds created_at & created_by
    # because ParamValue cannot be updated, we don't need updated_at
    created_at: datetime = DateTimeField(
        editable=False, db_default=models.functions.Now(), db_index=True
    )
    """Time of creation of record."""
    created_by: User = ForeignKey(
        User, PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""
    hash: str = CharField(max_length=HASH_LENGTH, null=True, db_index=True)

    class Meta:
        constraints = [
            # For simple types, use direct value comparison
            models.UniqueConstraint(
                fields=["param", "value"],
                name="unique_simple_param_value",
                condition=Q(hash__isnull=True),
            ),
            # For complex types (dictionaries), use hash
            models.UniqueConstraint(
                fields=["param", "hash"],
                name="unique_complex_param_value",
                condition=Q(hash__isnull=False),
            ),
        ]

    @classmethod
    def get_or_create(cls, param, value):
        # Simple types: int, float, str, bool
        if isinstance(value, (int, float, str, bool)):
            try:
                return cls.objects.create(param=param, value=value, hash=None), False
            except IntegrityError:
                return cls.objects.get(param=param, value=value), True

        # Complex types: dict, list
        else:
            hash = hash_dict(value)
            try:
                return cls.objects.create(param=param, value=value, hash=hash), False
            except IntegrityError:
                return cls.objects.get(param=param, hash=hash), True
