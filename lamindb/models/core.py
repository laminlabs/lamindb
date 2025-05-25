from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    overload,
)

from django.db import models

from lamindb.base.fields import (
    CharField,
)

from ..base.ids import base62_12
from .run import TracksRun, TracksUpdates
from .sqlrecord import SQLRecord

if TYPE_CHECKING:
    from pathlib import Path

    from upath import UPath

    from .artifact import Artifact


class Storage(SQLRecord, TracksRun, TracksUpdates):
    """Storage locations of artifacts such as S3 buckets or local directories.

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

    class Meta(SQLRecord.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "root"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        editable=False, unique=True, max_length=12, default=base62_12, db_index=True
    )
    """Universal id, valid across DB instances."""
    root: str = CharField(db_index=True, unique=True)
    """Root path of storage (cloud or local path)."""
    description: str | None = CharField(db_index=True, null=True)
    """A description of what the storage location is used for (optional)."""
    type: str = CharField(max_length=30, db_index=True)
    """Can be "local" vs. "s3" vs. "gs"."""
    region: str | None = CharField(max_length=64, db_index=True, null=True)
    """Cloud storage region, if applicable."""
    instance_uid: str | None = CharField(max_length=12, db_index=True, null=True)
    """Instance that manages this storage location."""
    artifacts: Artifact
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
