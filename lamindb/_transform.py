from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from lamin_utils import logger
from lamindb_setup.core._docs import doc_args

from lamindb.models import Run, Transform

from ._parents import _view_parents
from ._run import delete_run_artifacts
from .core._settings import settings
from .core.exceptions import InconsistentKey
from .core.versioning import message_update_key_in_version_family, process_revises

if TYPE_CHECKING:
    from lamindb.base.types import TransformType


def __init__(transform: Transform, *args, **kwargs):
    if len(args) == len(transform._meta.concrete_fields):
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    key: str | None = kwargs.pop("key") if "key" in kwargs else None
    description: str | None = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    revises: Transform | None = kwargs.pop("revises") if "revises" in kwargs else None
    version: str | None = kwargs.pop("version") if "version" in kwargs else None
    type: TransformType | None = kwargs.pop("type") if "type" in kwargs else "pipeline"
    reference: str | None = kwargs.pop("reference") if "reference" in kwargs else None
    reference_type: str | None = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    using_key = (
        kwargs.pop("using_key") if "using_key" in kwargs else settings._using_key
    )
    if "name" in kwargs:
        if key is None:
            key = kwargs.pop("name")
            warnings.warn(
                f"`name` will be removed soon, please pass '{key}' to `key` instead",
                FutureWarning,
                stacklevel=2,
            )
        else:
            # description wasn't exist, so no check necessary
            description = kwargs.pop("name")
            warnings.warn(
                f"`name` will be removed soon, please pass '{description}' to `description` instead",
                FutureWarning,
                stacklevel=2,
            )
    # below is internal use that we'll hopefully be able to eliminate
    uid: str | None = kwargs.pop("uid") if "uid" in kwargs else None
    if not len(kwargs) == 0:
        raise ValueError(
            "Only key, description, version, type, revises, reference, "
            f"reference_type can be passed, but you passed: {kwargs}"
        )
    if revises is None:
        # need to check uid before checking key
        if uid is not None:
            revises = (
                Transform.objects.using(using_key)
                .filter(uid__startswith=uid[:-4], is_latest=True)
                .order_by("-created_at")
                .first()
            )
        elif key is not None:
            candidate_for_revises = (
                Transform.objects.using(using_key)
                .filter(key=key, is_latest=True)
                .order_by("-created_at")
                .first()
            )
            if candidate_for_revises is not None:
                revises = candidate_for_revises
                if candidate_for_revises.source_code is None:
                    # no source code was yet saved, return the same transform
                    logger.important(
                        "no source code was yet saved, returning existing transform with same key"
                    )
                    uid = revises.uid
    if revises is not None and uid is not None and uid == revises.uid:
        from ._record import init_self_from_db, update_attributes

        if revises.key != key:
            logger.warning("ignoring inconsistent key")
        init_self_from_db(transform, revises)
        update_attributes(transform, {"description": description})
        return None
    if revises is not None and key is not None and revises.key != key:
        note = message_update_key_in_version_family(
            suid=revises.stem_uid,
            existing_key=revises.key,
            new_key=key,
            registry="Transform",
        )
        raise InconsistentKey(
            f"`key` is '{key}', but `revises.key` is '{revises.key}'\n\nEither do *not* pass `key`.\n\n{note}"
        )
    new_uid, version, key, description, revises = process_revises(
        revises, version, key, description, Transform
    )
    # this is only because the user-facing constructor allows passing a uid
    # most others don't
    if uid is None:
        has_consciously_provided_uid = False
        uid = new_uid
    else:
        has_consciously_provided_uid = True
    super(Transform, transform).__init__(
        uid=uid,
        description=description,
        key=key,
        type=type,
        version=version,
        reference=reference,
        reference_type=reference_type,
        _has_consciously_provided_uid=has_consciously_provided_uid,
        revises=revises,
    )


def delete(self) -> None:
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
        field="key",
        with_children=with_successors,
        distance=distance,
        attr_name="predecessors",
    )


Transform.__init__ = __init__
Transform.delete = delete
Transform.latest_run = latest_run
Transform.view_lineage = view_lineage
