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
