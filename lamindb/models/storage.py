from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    overload,
)

from django.db import models
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core._hub_core import (
    delete_storage_record,
    get_storage_records_for_instance,
)
from lamindb_setup.core._settings_storage import StorageSettings, init_storage
from lamindb_setup.core.upath import check_storage_is_empty, create_path

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
    """Storage locations of artifacts such as folders and S3 buckets.

    A storage location is either a folder (local or in the cloud) or
    an entire S3/GCP bucket.

    A LaminDB instance can manage and link multiple storage locations. But any
    storage location is managed by *at most one* LaminDB instance.

    .. dropdown:: Managed vs. referenced storage locations

        The LaminDB instance can write artifacts in managed storage
        locations but only read artifacts in referenced storage locations.

        The :attr:`~lamindb.Storage.instance_uid` field defines the managing LaminDB instance of a
        storage location. Some storage locations may not be managed by any LaminDB
        instance, in which case the `instance_uid` is `None`.

        When you delete a LaminDB instance, you'll be warned about data in managed
        storage locations while data in referenced storage locations is ignored.

    See Also:
        :attr:`lamindb.core.Settings.storage`
            Current default storage location of your compute session for writing artifacts.
        :attr:`~lamindb.setup.core.StorageSettings`
            Storage settings.

    Examples:

        When you create a LaminDB instance, you configure its default storage location via `--storage`::

            lamin init --storage ./myfolder  # or "s3://my-bucket/myfolder" or "gs://my-bucket/myfolder"

        View the current default storage location in your compute session for writing artifacts::

            import lamindb as ln

            ln.settings.storage

        Switch to another default storage location for writing artifacts::

            ln.settings.storage = "./myfolder2"  # or "s3://my-bucket/my-folder2" or "gs://my-bucket/my-folder2"

        View all storage locations used in your LaminDB instance::

            ln.Storage.df()

        Create a new storage location::

            ln.Storage(root="./myfolder3").save()

    Notes:

        .. dropdown:: How do I manage access to a storage location?

            You can low-level manage access through AWS policies that you attach to your S3 bucket
            or leverage LaminHub's fine-grained access management.

            :doc:`docs:access` explains both approaches.

        .. dropdown:: What is the `.lamindb/` directory inside a storage location?

            It stores all artifacts that are ingested through `lamindb`, indexed by the artifact `uid`.
            This means you don't have to worry about renaming or moving files, as this all happens on the database level.

            Existing artifacts are typically stored in hierarchical structures with semantic folder names.
            Instead of copying such artifacts into `.lamindb/` upon calls of `Artifact("legacy_path").save()`,
            LaminDB registers them with the semantic `key` representing the relative path within the storage location.
            These artifacts are marked with `artifact._key_is_virtual = False` and treated correspondingly.

            There is only a single `.lamindb/` directory per storage location.

        .. dropdown:: What should I do if I want to bulk migrate all artifacts to another storage?

            Currently, you can only achieve this manually and you should be careful with it.

            1. Copy or move artifacts into the desired new storage location
            2. Adapt the corresponding record in the {class}`~lamindb.Storage` registry by setting the `root` field to the new location
            3. If your LaminDB storage location is managed through the hub, you also need to update the storage record on the hub -- contact support

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
        if len(args) == len(self._meta.concrete_fields):
            super().__init__(*args)
            return None
        storage_record = Storage.filter(root=kwargs["root"]).one_or_none()
        if storage_record is not None:
            from .sqlrecord import init_self_from_db

            init_self_from_db(self, storage_record)
            return None
        if "_skip_preparation" in kwargs:
            skip_preparation = kwargs.pop("_skip_preparation")
            if skip_preparation is True:
                super().__init__(*args, **kwargs)
                return None
            kwargs.pop("_skip_preparation")

        # instance_id won't take effect if
        # - there is no write access
        # - the storage location is already managed by another instance
        ssettings, _ = init_storage(
            kwargs["root"],
            instance_id=setup_settings.instance._id,
            prevent_register_hub=not setup_settings.instance.is_on_hub,
        )
        assert kwargs["root"] == ssettings.root_as_str  # noqa: S101
        if "instance_uid" in kwargs:
            assert kwargs["instance_uid"] == ssettings.instance_uid  # noqa: S101
        else:
            kwargs["instance_uid"] = ssettings.instance_uid
        if ssettings._uid is not None:  # need private attribute here
            kwargs["uid"] = ssettings._uid
        if "type" not in kwargs:
            kwargs["type"] = ssettings.type
        else:
            assert kwargs["type"] == ssettings.type  # noqa: S101
        if "region" in kwargs:
            assert kwargs["region"] == ssettings.region  # noqa: S101
        else:
            kwargs["region"] = ssettings.region

        if ssettings.instance_uid is not None:
            hub_message = ""
            if setup_settings.instance.is_on_hub:
                instance_owner = setup_settings.instance.owner
                hub_message = f", see: https://lamin.ai/{instance_owner}/infrastructure"
            managed_message = (
                "managed" if ssettings.instance_uid else "referenced (read-only)"
            )
            logger.important(f"created {managed_message} storage location{hub_message}")
        super().__init__(**kwargs)

    @property
    def path(self) -> Path | UPath:
        """Path.

        Uses the `.root` field and converts it into a `Path` or `UPath`.
        """
        access_token = self._access_token if hasattr(self, "_access_token") else None
        return create_path(self.root, access_token=access_token)

    def delete(self) -> None:
        """Delete the storage location.

        This errors in case the storage location is not empty.
        """
        from .. import settings

        assert not self.artifacts.exists(), "Cannot delete storage holding artifacts."  # noqa: S101
        check_storage_is_empty(self.path)
        assert settings.storage.root_as_str != self.root, (  # noqa: S101
            "Cannot delete the current storage location, switch to another."
        )
        if setup_settings.user.handle != "anonymous":  # only attempt if authenticated
            storage_records = get_storage_records_for_instance(
                setup_settings.instance._id
            )
            for storage_record in storage_records:
                if storage_record["lnid"] == self.uid:
                    assert storage_record["is_default"] in {False, None}, (  # noqa: S101
                        "Cannot delete default storage of instance."
                    )
                    delete_storage_record(storage_record)
        ssettings = StorageSettings(self.root)
        if ssettings._mark_storage_root.exists():
            ssettings._mark_storage_root.unlink(
                missing_ok=True  # this is totally weird, but needed on Py3.11
            )
        super().delete()
