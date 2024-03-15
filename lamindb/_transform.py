from typing import Optional

from lnschema_core.models import Artifact, Run, Transform
from lnschema_core.types import TransformType

from ._run import delete_run_artifacts
from .core.versioning import get_uid_from_old_version, init_uid


def __init__(transform: Transform, *args, **kwargs):
    if len(args) == len(transform._meta.concrete_fields):
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    key: Optional[str] = kwargs.pop("key") if "key" in kwargs else None
    is_new_version_of: Optional[Transform] = (
        kwargs.pop("is_new_version_of") if "is_new_version_of" in kwargs else None
    )
    (kwargs.pop("initial_version_id") if "initial_version_id" in kwargs else None)
    version: Optional[str] = kwargs.pop("version") if "version" in kwargs else None
    type: Optional[TransformType] = (
        kwargs.pop("type") if "type" in kwargs else TransformType.pipeline
    )
    reference: Optional[str] = (
        kwargs.pop("reference") if "reference" in kwargs else None
    )
    reference_type: Optional[str] = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    # below is internal use that we'll hopefully be able to eliminate
    uid: Optional[str] = kwargs.pop("uid") if "uid" in kwargs else None
    if not len(kwargs) == 0:
        raise ValueError(
            "Only name, key, version, type, is_new_version_of, reference, "
            f"reference_type can be passed, but you passed: {kwargs}"
        )
    if is_new_version_of is None:
        new_uid = init_uid(version=version, n_full_id=Transform._len_full_uid)
    else:
        if not isinstance(is_new_version_of, Transform):
            raise TypeError("is_new_version_of has to be of type ln.Transform")
        new_uid, version = get_uid_from_old_version(is_new_version_of, version)
        if name is None:
            name = is_new_version_of.name

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
    # set latest_report to None, it's tracked through the latest run
    latest_report = None
    if self.latest_report is not None:
        latest_report = self.latest_report
        self.latest_report = None
    source_code = None
    if self.source_code is not None:
        source_code = self.source_code
        self.source_code = None
    if latest_report is not None or source_code is not None:
        self.save()
    if source_code is not None:
        source_code.delete(permanent=True)
    # query all runs and delete their artifacts
    runs = Run.filter(transform=self)
    for run in runs:
        delete_run_artifacts(run)
    # at this point, all artifacts have been taken care of
    # we can now leverage CASCADE delete
    super(Transform, self).delete()


@property  # type: ignore
def latest_run(self) -> Run:
    return self.runs.order_by("-started_at").first()


Transform.__init__ = __init__
Transform.delete = delete
Transform.latest_run = latest_run
