from django.db import models

DEFAULT_CREATED_BY_UID = "0" * 8
DEFAULT_BRANCH_CODE = 1
DEFAULT_RUN_UID = "0" * 16


class WriteLogTableState(models.Model):
    """A list of tables for which we're recording write logs.

    This table serves two purposes: it allows us to store the
    name of the table for a given write log record with less state.

    NOTE: This table only makes sense for a single LaminDB instance,
    since the ID of a given table may be different on different instances.
    You should join this table into WriteLog to determine which table
    a given write log record applies to, rather than using the IDs directly.
    """

    # Using a smallserial here because it's incredibly unlikely that
    # we'll have more than 32k tables
    id = models.SmallAutoField(primary_key=True)
    table_name = models.CharField(max_length=255)
    backfilled = models.BooleanField()


class WriteLogMigrationState(models.Model):
    """A summary of the state of Django's migrations when a write log record was recorded.

    When a write log record is recorded, we need to record the state of the data migrations
    for each Django app so that we know if we need to migrate the record from the version
    it had when it was written to the version that exists in the target during replay.
    """

    # Using a smallserial here because it's incredibly unlikely that
    # we'll have more than 32k schema migrations.
    id = models.SmallAutoField(primary_key=True)
    migration_state_id = models.JSONField()


class WriteLog(models.Model):
    """Stores the write log for LaminDB tables."""

    id = models.BigAutoField(primary_key=True)
    migration_state = models.ForeignKey(
        WriteLogMigrationState, on_delete=models.PROTECT
    )
    table = models.ForeignKey(WriteLogTableState, on_delete=models.PROTECT)
    uid = models.CharField(max_length=18, editable=False, db_index=True, unique=True)
    # While all normal tables will have a space ID, many-to-many tables won't.
    space_uid = models.CharField(max_length=12, null=True)
    created_by_uid = models.CharField(max_length=8, default=DEFAULT_CREATED_BY_UID)
    branch_code = models.IntegerField(default=DEFAULT_BRANCH_CODE)
    run_uid = models.CharField(max_length=16, default=DEFAULT_RUN_UID)
    # Many-to-many tables don't have row UIDs, so this needs to be nullable.
    record_uid = models.JSONField(null=True)
    record_data = models.JSONField(null=True)
    event_type = models.PositiveSmallIntegerField()
    created_at = models.DateTimeField()

    class Meta:
        verbose_name = "Write Log"
        verbose_name_plural = "Write Logs"


class WriteLogLock(models.Model):
    """A lock table for write log triggers.

    Write log triggers for a given table record all changes to
    that table. We need to make sure that when we're replaying
    the write log from another node, we don't record the replayed
    changes as additional write log records. We'll accomplish that in a
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
