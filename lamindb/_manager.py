from typing import Optional

from django.db import models


class Manager(models.Manager):
    """Extension of Django Manager.

    See Also:

        `django Manager <https://docs.djangoproject.com/en/4.2/topics/db/managers/>`__

    Examples:

        >>> ln.save(ln.Tag.from_values(["Tag1", "Tag2", "Tag3"], field="name"))
        >>> tags = ln.Tag.select(name__icontains = "tag").all()
        >>> ln.Project(name="Project1").save()
        >>> project = ln.Project.select(name="Project1").one()
        >>> project.tags.set(tags)
        >>> manager = project.tags
    """

    def list(self, field: Optional[str] = None):
        """Populate a list with the results.

        Examples:
            >>> ln.save(ln.Tag.from_values(["Tag1", "Tag2", "Tag3"], field="name"))
            >>> tags = ln.Tag.select(name__icontains = "tag").all()
            >>> ln.Project(name="Project1").save()
            >>> project = ln.Project.select(name="Project1").one()
            >>> project.tags.set(tags)
            >>> project.tags.list()
            [Tag(id=sFMcPepC, name=Tag1, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            Tag(id=2SscQvsM, name=Tag2, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse), # noqa
            Tag(id=lecV87vi, name=Tag3, updated_at=2023-07-19 19:45:17, created_by_id=DzTjkKse)] # noqa
            >>> project.tags.list("name")
            ['Tag1', 'Tag2', 'Tag3']
        """
        if field is None:
            return [item for item in self.all()]
        else:
            return [item for item in self.values_list(field, flat=True)]


setattr(models.Manager, "list", Manager.list)
