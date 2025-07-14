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
from lamindb_setup.core._settings_storage import (
    StorageSettings,
    get_storage_type,
    init_storage,
)
from lamindb_setup.core.upath import check_storage_is_empty, create_path

from lamindb.base.fields import (
    CharField,
)

from ..base.ids import base62_12
from .run import TracksRun, TracksUpdates
from .sqlrecord import SQLRecord

if TYPE_CHECKING:
    from pathlib import Path

    from lamindb_setup.types import StorageType
    from upath import UPath

    from .artifact import Artifact


class Storage(SQLRecord, TracksRun, TracksUpdates):
    """Storage locations of artifacts such as local directories or S3 buckets.

    A storage location is either a directory (local or a folder in the cloud) or
    an entire S3/GCP bucket.
    A LaminDB instance can manage and read from multiple storage locations. But any
    storage location is managed by *at most one* LaminDB instance.

    .. dropdown:: Managed vs. read-only storage locations

        A LaminDB instance can only write artifacts to its managed storage
        locations.

        The :attr:`~lamindb.Storage.instance_uid` field defines the managing LaminDB instance of a storage location.
        You can access the `instance_uid` of your current instance through `ln.setup.settings.instance_uid`.

        Here is an example (`source <https://lamin.ai/laminlabs/lamindata/transform/dPco79GYgzag0000>`__).

        .. image:: https://lamin-site-assets.s3.amazonaws.com/.lamindb/eHDmIOAxLEoqZ2oK0000.png
           :width: 400px

        Some public storage locations are not be managed by any LaminDB instance: their `instance_uid` is `None`.

    .. dropdown:: Managing access to storage locations across instances

        You can manage access through AWS policies that you attach to your S3 bucket
        or leverage LaminHub's fine-grained access management.

        Head over to `https://lamin.ai/{account}/infrastructure`.
        By clicking the green button that says "Connect S3 bucket", you enable Lamin to issue federated S3 tokens
        for a bucket so that your collaborators can access data based on their permissions in LaminHub.
        :doc:`docs:access` has more details.

        .. image:: https://lamin-site-assets.s3.amazonaws.com/.lamindb/ze8hkgVxVptSSZEU0000.png
           :width: 800px

        If you don't want to store data in the cloud, you can use local storage locations: :doc:`faq/keep-artifacts-local`.

    Args:
        root: `str` The root path of the storage location, e.g., `"./mydir"`, `"s3://my-bucket"`, `"s3://my-bucket/myfolder"`, `"gs://my-bucket/myfolder"`, `"/nfs/shared/datasets/genomics"`, `"/weka/shared/models/"`, ...
        description: `str | None = None` An optional description.
        host: `str | None = None` For local storage locations, pass a globally unique host identifier, e.g. `"my-institute-cluster-1"`, `"my-server-abcd"`, ...

    See Also:
        :attr:`lamindb.core.Settings.storage`
            Current default storage location of your compute session for writing artifacts.
        :attr:`~lamindb.setup.core.StorageSettings`
            Storage settings.
        :doc:`faq/keep-artifacts-local`
            Avoid storing artifacts in the cloud, but keep them on local infrastructure.

    Examples:

        When you create a LaminDB instance, you configure its default storage location via `--storage`::

            lamin init --storage ./mydatadir  # or "s3://my-bucket/myfolder", "gs://my-bucket/myfolder", ...

        View the current default storage location for writing artifacts::

            import lamindb as ln

            ln.settings.storage

        Create a new cloud storage location::

            ln.Storage(root="s3://our-bucket/our-folder").save()

        Create a new local storage location::

            ln.Storage(root="/dir/our-shared-dir", host="our-server-123").save()

        Switch to another storage location::

            ln.settings.storage = "/dir/our-shared-dir"  # or "s3://our-bucket/our-folder", "gs://our-bucket/our-folder", ...

        If you're operating in `keep-artifacts-local` mode (:doc:`faq/keep-artifacts-local`), you can switch among additional local storage locations::

            ln.Storage(root="/dir/our-other-shared-dir", host="our-server-123").save()  # create
            ln.settings.local_storage = "/dir/our-other-shared-dir"  # switch

        View all storage locations used in your LaminDB instance::

            ln.Storage.df()

    Notes:

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
    type: StorageType = CharField(max_length=30, db_index=True)
    """Can be "local" vs. "s3" vs. "gs". Is auto-detected from the format of the `root` path."""
    region: str | None = CharField(max_length=64, db_index=True, null=True)
    """Storage region for cloud storage locations. Host identifier for local storage locations."""
    instance_uid: str | None = CharField(max_length=12, db_index=True, null=True)
    """Instance that manages this storage location."""
    artifacts: Artifact
    """Artifacts contained in this storage location."""

    @overload
    def __init__(
        self,
        root: str,
        *,
        description: str | None = None,
        host: str | None = None,
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
        if args:
            assert len(args) == 1, (  # noqa: S101
                "Storage can only be initialized with a single positional argument, the root path."
            )
            kwargs["root"] = args[0]
        if "host" in kwargs:
            if "type" in kwargs:
                assert kwargs["type"] == "local", (  # noqa: S101
                    "type needs to be 'local' if host is set"
                )
            else:
                kwargs["type"] = "local"
            assert get_storage_type(kwargs["root"]) == "local", (  # noqa: S101
                "root must be a local path if host is set"
            )
            assert "region" not in kwargs, "region must not be set if host is set"  # noqa: S101
            kwargs["region"] = kwargs.pop("host")
            storage_record = Storage.filter(
                root=kwargs["root"], region=kwargs["region"]
            ).one_or_none()
        else:
            storage_record = Storage.filter(root=kwargs["root"]).one_or_none()
        if storage_record is not None:
            from .sqlrecord import init_self_from_db

            init_self_from_db(self, storage_record)
            return None

        skip_preparation = kwargs.pop("_skip_preparation", False)
        if skip_preparation:
            super().__init__(*args, **kwargs)
            return None

        # instance_id won't take effect if
        # - there is no write access
        # - the storage location is already managed by another instance
        ssettings, _ = init_storage(
            kwargs["root"],
            instance_id=setup_settings.instance._id,
            instance_slug=setup_settings.instance.slug,
            register_hub=setup_settings.instance.is_on_hub,
            prevent_register_hub=not setup_settings.instance.is_on_hub,
            region=kwargs.get("region", None),  # host was renamed to region already
        )
        # ssettings performed validation and normalization of the root path
        kwargs["root"] = ssettings.root_as_str  # noqa: S101
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

        is_managed_by_current_instance = (
            ssettings.instance_uid == setup_settings.instance.uid
        )
        if ssettings.instance_uid is not None and not is_managed_by_current_instance:
            is_managed_by_instance = (
                f", is managed by instance with uid {ssettings.instance_uid}"
            )
        else:
            is_managed_by_instance = ""
        hub_message = ""
        if setup_settings.instance.is_on_hub and is_managed_by_current_instance:
            instance_owner = setup_settings.instance.owner
            hub_message = f", see: https://lamin.ai/{instance_owner}/infrastructure"
        managed_message = (
            "created managed"
            if is_managed_by_current_instance
            else "referenced read-only"
        )
        logger.important(
            f"{managed_message} storage location at {kwargs['root']}{is_managed_by_instance}{hub_message}"
        )
        super().__init__(**kwargs)

    @property
    def host(self) -> str | None:
        """Host identifier for local storage locations.

        Is `None` for locations with `type != "local"`.

        A globally unique user-defined host identifier (cluster, server, laptop, etc.).
        """
        if self.type != "local":
            return None
        return self.region

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
                # only query those storage records on the hub that are managed by the current instance
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
