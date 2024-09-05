from __future__ import annotations

from typing import TYPE_CHECKING

from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import Run, Transform

from ._parents import _view_parents
from ._run import delete_run_artifacts
from .core.versioning import message_update_key_in_version_family, process_revises

if TYPE_CHECKING:
    from lnschema_core.types import TransformType


def __init__(transform: Transform, *args, **kwargs):
    if len(args) == len(transform._meta.concrete_fields):
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    name: str | None = kwargs.pop("name") if "name" in kwargs else None
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
    revises: Transform | None = kwargs.pop("revises") if "revises" in kwargs else None
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    type: TransformType | None = kwargs.pop("type") if "type" in kwargs else "pipeline"
    reference: str | None = kwargs.pop("reference") if "reference" in kwargs else None
    reference_type: str | None = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    if "is_new_version_of" in kwargs:
        logger.warning("`is_new_version_of` will be removed soon, please use `revises`")
        revises = kwargs.pop("is_new_version_of")
    # below is internal use that we'll hopefully be able to eliminate
    uid: str | None = kwargs.pop("uid") if "uid" in kwargs else None
    if not len(kwargs) == 0:
        raise ValueError(
            "Only name, key, version, type, revises, reference, "
            f"reference_type can be passed, but you passed: {kwargs}"
        )
    if revises is None:
        if key is not None:
            revises = Transform.filter(key=key).order_by("-created_at").first()
        elif uid is not None and not uid.endswith("0000"):
            revises = Transform.filter(uid__startswith=uid[:-4]).one_or_none()
    if revises is not None and key is not None and revises.key != key:
        note = message_update_key_in_version_family(
            suid=revises.stem_uid,
            existing_key=revises.key,
            new_key=key,
            registry="Artifact",
        )
        raise ValueError(
            f"`key` is {key}, but `revises.key` is '{revises.key}'\n\n Either do *not* pass `key`.\n\n{note}"
        )
    new_uid, version, name, revises = process_revises(revises, version, name, Transform)
    # this is only because the user-facing constructor allows passing a uid
    # most others don't
    if uid is None:
        has_consciously_provided_uid = False
        uid = new_uid
    else:
        has_consciously_provided_uid = True
    super(Transform, transform).__init__(
        uid=uid,
        name=name,
        key=key,
        type=type,
        version=version,
        reference=reference,
        reference_type=reference_type,
        _has_consciously_provided_uid=has_consciously_provided_uid,
        revises=revises,
    )


def delete(self) -> None:
    _source_code_artifact = None
    if self._source_code_artifact is not None:
        _source_code_artifact = self._source_code_artifact
        self._source_code_artifact = None
        self.save()
    if _source_code_artifact is not None:
        _source_code_artifact.delete(permanent=True)
    # query all runs and delete their artifacts
    runs = Run.filter(transform=self)
    for run in runs:
        delete_run_artifacts(run)
    # at this point, all artifacts have been taken care of
    # we can now leverage CASCADE delete
    super(Transform, self).delete()


@property  # type: ignore
@doc_args(Transform.latest_run.__doc__)
def latest_run(self) -> Run:
    """{}"""  # noqa: D415
    return self.runs.order_by("-started_at").first()


def view_lineage(self, with_successors: bool = False, distance: int = 5):
    return _view_parents(
        record=self,
        field="name",
        with_children=with_successors,
        distance=distance,
        attr_name="predecessors",
    )


Transform.__init__ = __init__
Transform.delete = delete
Transform.latest_run = latest_run
Transform.view_lineage = view_lineage
