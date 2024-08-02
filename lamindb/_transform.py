from __future__ import annotations

from typing import TYPE_CHECKING

from lamindb_setup.core._docs import doc_args
from lnschema_core.models import Run, Transform

from ._parents import _view_parents
from ._run import delete_run_artifacts
from .core.versioning import process_is_new_version_of

if TYPE_CHECKING:
    from lnschema_core.types import TransformType


def __init__(transform: Transform, *args, **kwargs):
    if len(args) == len(transform._meta.concrete_fields):
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    name: str | None = kwargs.pop("name") if "name" in kwargs else None
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
    is_new_version_of: Transform | None = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    (kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None)
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    type: TransformType | None = kwargs.pop("type") if "type" in kwargs else "pipeline"
    reference: str | None = kwargs.pop("reference") if "reference" in kwargs else None
    reference_type: str | None = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    # below is internal use that we'll hopefully be able to eliminate
    uid: str | None = kwargs.pop("uid") if "uid" in kwargs else None
    if not len(kwargs) == 0:
        raise ValueError(
            "Only name, key, version, type, is_new_version_of, reference, "
            f"reference_type can be passed, but you passed: {kwargs}"
        )
    new_uid, version, name = process_is_new_version_of(
        is_new_version_of, version, name, Transform
    )
    # this is only because the user-facing constructor allows passing an id
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
        record=self, field=None, with_children=with_successors, distance=distance
    )


Transform.__init__ = __init__
Transform.delete = delete
Transform.latest_run = latest_run
Transform.view_lineage = view_lineage
