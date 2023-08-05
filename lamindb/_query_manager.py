from typing import Optional

from django.db import models


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


setattr(models.Manager, "list", QueryManager.list)
setattr(models.Manager, "df", QueryManager.df)
