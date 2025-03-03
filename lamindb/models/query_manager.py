from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

from django.db import models
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args

from lamindb.models import Record

from ..core._settings import settings

if TYPE_CHECKING:
    from lamindb.base.types import StrField


class QueryManager(models.Manager):
    """Manage queries through fields.

    See Also:

        :class:`lamindb.core.QuerySet`
        `django Manager <https://docs.djangoproject.com/en/4.2/topics/db/managers/>`__

    Examples:

        >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))  # noqa
        >>> labels = ln.ULabel.filter(name__icontains = "label").all()
        >>> ln.ULabel(name="ULabel1").save()
        >>> label = ln.ULabel.get(name="ULabel1")
        >>> label.parents.set(labels)
        >>> manager = label.parents
        >>> manager.df()
    """

    def _track_run_input_manager(self):
        if hasattr(self, "source_field_name") and hasattr(self, "target_field_name"):
            if (
                self.source_field_name == "collection"
                and self.target_field_name == "artifact"
            ):
                from lamindb.core._context import context
                from lamindb.models.artifact import (
                    WARNING_RUN_TRANSFORM,
                    _track_run_input,
                )

                if (
                    context.run is None
                    and not settings.creation.artifact_silence_missing_run_warning
                ):
                    logger.warning(WARNING_RUN_TRANSFORM)
                _track_run_input(self.instance)

    def list(self, field: str | None = None):
        """Populate a list with the results.

        Examples:
            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))
            >>> labels = ln.ULabel.filter(name__icontains="label").all()
            >>> ln.ULabel(name="ULabel1").save()
            >>> label = ln.ULabel.get(name="ULabel1")
            >>> label.parents.set(labels)
            >>> label.parents.list()
            >>> label.parents.list("name")
            ['ULabel1', 'ULabel2', 'ULabel3']
        """
        self._track_run_input_manager()
        if field is None:
            return list(self.all())
        else:
            return list(self.values_list(field, flat=True))

    def df(self, **kwargs):
        """Convert to DataFrame.

        For `**kwargs`, see :meth:`lamindb.core.QuerySet.df`.
        """
        return self.all().df(**kwargs)

    def all(self):
        """Return QuerySet of all.

        For `**kwargs`, see :meth:`lamindb.core.QuerySet.df`.
        """
        self._track_run_input_manager()
        return self._all_base_class()

    @doc_args(Record.search.__doc__)
    def search(self, string: str, **kwargs):
        """{}"""  # noqa: D415
        from .record import _search

        return _search(cls=self.all(), string=string, **kwargs)

    @doc_args(Record.lookup.__doc__)
    def lookup(self, field: StrField | None = None, **kwargs) -> NamedTuple:
        """{}"""  # noqa: D415
        from .record import _lookup

        return _lookup(cls=self.all(), field=field, **kwargs)


models.Manager.list = QueryManager.list
models.Manager.df = QueryManager.df
models.Manager.search = QueryManager.search
models.Manager.lookup = QueryManager.lookup
models.Manager._track_run_input_manager = QueryManager._track_run_input_manager
# the two lines below would be easy if we could actually inherit; like this,
# they're suboptimal
models.Manager._all_base_class = models.Manager.all
models.Manager.all = QueryManager.all
