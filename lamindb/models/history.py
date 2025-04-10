from django.db import models


class History(models.Model):
    """Stores the edit history for all LaminDB tables."""

    seqno = models.AutoField(primary_key=True)
    migration_history_id = models.BigIntegerField()
    id = models.UUIDField()
    table_name = models.CharField(max_length=100)
    record_id = models.BigIntegerField()
    record_data = models.JSONField(null=True)
    event_type = models.PositiveSmallIntegerField()
    created_at = models.DateTimeField(auto_now_add=True)

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
