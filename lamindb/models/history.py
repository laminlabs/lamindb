from datetime import datetime
from uuid import UUID

# from .record import User --> only needed for ForeignKey
# from django.contrib.contenttypes.models import ContentType --> only needed for ForeignKey
from django.db.models import (
    BigIntegerField,
    DateTimeField,
    JSONField,
    SmallIntegerField,
    UUIDField,
    functions,
)

from .record import BasicRecord


class History(BasicRecord):
    id: UUID = UUIDField(
        primary_key=True
    )  # snowflake ID would be ideal because then we'd have some ordering; otherwise timsstamped uuid also ok
    table_id: int = SmallIntegerField()  # ForeignKey into ContentType but unprotected so that tables can be deleted; so not using ForeignKey
    operation: int = (
        SmallIntegerField()
    )  # code INSERT, DELETE, UPDATE and other things we might want
    record_id: int = BigIntegerField()
    values: dict = JSONField()
    created_at: datetime = DateTimeField(
        db_default=functions.Now(),
    )  # not strictly necessary if id is timestamped
    created_by: int = SmallIntegerField()  # ForeignKey into User but unprotected so that we don't pay the price for the constraint; no need for integrity as table is immutable
