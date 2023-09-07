from typing import Optional

from django.db import models

from .dev._feature_manager import get_feature_set_by_slot


class QueryManager(models.Manager):
    """Manage queries through fields.

    See Also:

        :class:`lamindb.dev.QuerySet`
        `django Manager <https://docs.djangoproject.com/en/4.2/topics/db/managers/>`__

    Examples:

        >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))  # noqa
        >>> labels = ln.ULabel.filter(name__icontains = "label").all()
        >>> ln.ULabel(name="ULabel1").save()
        >>> label = ln.ULabel.filter(name="ULabel1").one()
        >>> label.parents.set(labels)
        >>> manager = label.parents
        >>> manager.df()
    """

    def list(self, field: Optional[str] = None):
        """Populate a list with the results.

        Examples:
            >>> ln.save(ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name"))
            >>> labels = ln.ULabel.filter(name__icontains = "label").all()
            >>> ln.ULabel(name="ULabel1").save()
            >>> label = ln.ULabel.filter(name="ULabel1").one()
            >>> label.parents.set(labels)
            >>> label.parents.list()
            [ULabel(id=sFMcPepC, name=ULabel1, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            ULabel(id=2SscQvsM, name=ULabel2, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            ULabel(id=lecV87vi, name=ULabel3, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse)] # noqa
            >>> label.parents.list("name")
            ['ULabel1', 'ULabel2', 'ULabel3']
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
