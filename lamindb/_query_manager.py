from typing import Dict, List, Optional, Union

from django.db import models
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import Data, Dataset, Feature, File, Registry

from ._query_set import QuerySet
from .dev._data import add_labels, get_labels
from .dev._feature_manager import get_feature_set_by_slot


class QueryManager(models.Manager):
    """Manage queries through fields.

    See Also:

        :class:`lamindb.dev.QuerySet`
        `django Manager <https://docs.djangoproject.com/en/4.2/topics/db/managers/>`__

    Examples:

        >>> ln.save(ln.Label.from_values(["Label1", "Label2", "Label3"], field="name"))
        >>> labels = ln.Label.filter(name__icontains = "label").all()
        >>> ln.Label(name="Label1").save()
        >>> label = ln.Label.filter(name="Label1").one()
        >>> label.parents.set(labels)
        >>> manager = label.parents
        >>> manager.df()
    """

    def list(self, field: Optional[str] = None):
        """Populate a list with the results.

        Examples:
            >>> ln.save(ln.Label.from_values(["Label1", "Label2", "Label3"], field="name"))
            >>> labels = ln.Label.filter(name__icontains = "label").all()
            >>> ln.Label(name="Label1").save()
            >>> label = ln.Label.filter(name="Label1").one()
            >>> label.parents.set(labels)
            >>> label.parents.list()
            [Label(id=sFMcPepC, name=Label1, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            Label(id=2SscQvsM, name=Label2, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            Label(id=lecV87vi, name=Label3, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse)] # noqa
            >>> label.parents.list("name")
            ['Label1', 'Label2', 'Label3']
        """
        if field is None:
            return [item for item in self.all()]
        else:
            return [item for item in self.values_list(field, flat=True)]

    def df(self, **kwargs):
        """Convert to DataFrame.

        For `**kwargs`, see :meth:`lamindb.dev.QuerySet.df`.
        """
        return self.all().df(**kwargs)

    @doc_args(Data.add_labels.__doc__)
    def add_by_feature(
        self,
        records: Union[Registry, List[Registry], QuerySet],
        feature: Feature,
    ) -> None:
        """{}"""
        if not isinstance(self.isinstance, (Dataset, File)):
            raise TypeError("Instance must be File or Dataset.")
        return add_labels(self=self.instance, records=records, feature=feature)

    @doc_args(Data.get_labels.__doc__)
    def get_by_feature(
        self,
        feature: Feature,
        mute: bool = False,
        flat_names: bool = False,
    ) -> Union[QuerySet, Dict[str, QuerySet], List]:
        """{}"""
        if not isinstance(self.isinstance, (Dataset, File)):
            raise TypeError("Instance must be File or Dataset.")
        return get_labels(
            self=self.instance, feature=feature, mute=mute, flat_names=flat_names
        )

    def __getitem__(self, item: str):
        try:
            source_field_name = self.source_field_name
            target_field_name = self.target_field_name

            if (
                source_field_name in {"file", "dataset"}
                and target_field_name == "feature_set"
            ):
                return get_feature_set_by_slot(host=self.instance).get(item)

        except Exception:  # pragma: no cover
            return


setattr(models.Manager, "list", QueryManager.list)
setattr(models.Manager, "df", QueryManager.df)
setattr(models.Manager, "__getitem__", QueryManager.__getitem__)
setattr(models.Manager, "add_by_feature", QueryManager.add_by_feature)
setattr(models.Manager, "get_by_feature", QueryManager.get_by_feature)
