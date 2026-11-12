from django.db import models


class ArtifactUIDReservation(models.Model):
    """Internal, durable allocation of a streamed artifact's managed path.

    Keep terminal reservations as tombstones: an abandoned object path must never
    be allocated to another writer.
    """

    uid = models.CharField(max_length=20, primary_key=True)
    storage = models.ForeignKey("lamindb.Storage", on_delete=models.SET_NULL, null=True)
    key = models.CharField(max_length=1024, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    state = models.CharField(max_length=20, default="reserved")

    class Meta:
        app_label = "lamindb"


def lock_artifact_uid(artifact, using):
    """Coordinate managed-path allocation inside the caller's transaction."""
    from lamindb.errors import IntegrityError

    from .storage import Storage

    Storage.objects.using(using).select_for_update().get(pk=artifact.storage_id)
    artifact.__class__.objects.using(using).select_for_update().filter(
        uid__startswith=artifact.stem_uid
    ).order_by("id").first()
    if (
        artifact._state.adding
        and getattr(artifact, "_streaming_reserved_uid", None) != artifact.uid
        and ArtifactUIDReservation.objects.using(using)
        .filter(uid=artifact.uid)
        .exists()
    ):
        raise IntegrityError(
            f"Artifact UID {artifact.uid} is reserved by a streaming writer. "
            "Reconstruct the artifact to allocate a new version."
        )
