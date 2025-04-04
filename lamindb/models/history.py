from uuid import UUID

from django.contrib.contenttypes.models import ContentType
from django.db.models import JSONField, SmallIntegerField, UUIDField

from .record import BasicRecord


class History(BasicRecord):
    id: UUID = UUIDField(
        primary_key=True
    )  # snowflake ID would be ideal because then we'd have some ordering; otherwise timsstamped uuid also ok
    table: ContentType = SmallIntegerField()  # ForeignKey into ContentType but unprotected so that tables can be deleted; so not using ForeignKey
    operation: int = (
        SmallIntegerField()
    )  # code INSERT, DELETE, UPDATE and other things we might want
    data: dict = JSONField()
