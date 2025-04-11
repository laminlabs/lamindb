from django.db import models


class HistoryTableState(models.Model):
    """A list of tables for which we're tracking history.

    This table serves two purposes: it allows us to store the
    name of the table for a given history record with less state.

    NOTE: This table only makes sense for a single LaminDB instance,
    since the ID of a given table may be different on different instances.
    You should join this table into History to determine which table
    a given history record applies to, rather than using the IDs directly.
    """

    # Using a smallserial here because it's incredibly unlikely that
    # we'll have more than 32k tables
    id = models.SmallAutoField(primary_key=True)
    table_name = models.CharField(max_length=255)
    backfilled = models.BooleanField()


class HistoryMigrationState(models.Model):
    """A summary of the state of Django's migrations when a history record was recorded.

    When a history record is recorded, we need to record the state of the data migrations
    for each Django app so that we know if we need to migrate the record from the version
    it had when it was written to the version that exists in the target during replay.
    """

    # Using a smallserial here because it's incredibly unlikely that
    # we'll have more than 32k schema migrations.
    id = models.SmallAutoField(primary_key=True)
    migration_history_id = models.JSONField()


class History(models.Model):
    """Stores the edit history for all LaminDB tables."""

    seqno = models.AutoField(primary_key=True)
    migration_history_id = models.ForeignKey(
        HistoryMigrationState, on_delete=models.PROTECT
    )
    table_id = models.ForeignKey(HistoryTableState, on_delete=models.PROTECT)
    id = models.UUIDField()
    # Many-to-many tables don't have row UIDs, so this needs to be nullable.
    record_uid = models.CharField(max_length=8, null=True)
    record_data = models.JSONField(null=True)
    event_type = models.PositiveSmallIntegerField()
    created_at = models.DateTimeField()

    class Meta:
        verbose_name = "Edit History"
        verbose_name_plural = "Edit Histories"


class HistoryLock(models.Model):
    """A lock table for history triggers.

    History triggers for a given table record all changes to
    that table. We need to make sure that when we're replaying
    history from another node, we don't record the replayed
    changes as additional history. We'll accomplish that in a
    database-agnostic way using this lock table.

    The table is presumed to have a single row, whose value
    is a boolean.
    """

    locked = models.BooleanField()

    def save(self, *args, **kwargs):
        # Ensure that there's only one object in the table.
        self.__class__.objects.exclude(id=self.pk).delete()
        super().save(*args, **kwargs)

    @classmethod
    def load(cls):
        if cls.objects.count() == 0:
            return cls.objects.create(locked=False)
        return cls.objects.first()

    def lock(self):
        self.locked = True
        self.save()

    def unlock(self):
        self.locked = False
        self.save()
