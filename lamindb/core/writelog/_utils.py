import re

from django.db import models

from lamindb.models.writelog import MigrationState, TableState


class DjangoMigration(models.Model):
    """This model class allows us to access the migrations table using normal Django syntax."""

    app = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    applied = models.DateTimeField()

    class Meta:
        managed = False  # Tell Django not to manage this table
        db_table = "django_migrations"  # Specify the actual table name


def update_write_log_table_state(tables: set[str]):
    existing_tables = set(
        TableState.objects.filter(table_name__in=tables).values_list(
            "table_name", flat=True
        )
    )

    table_states_to_add = [
        TableState(table_name=t, backfilled=False)
        for t in tables
        if t not in existing_tables
    ]

    TableState.objects.bulk_create(table_states_to_add)


def get_latest_migration_state() -> MigrationState | None:
    return MigrationState.objects.order_by("-id").first()


def update_migration_state():
    app_migrations = {}

    for app in DjangoMigration.objects.values_list("app", flat=True).distinct():
        migrations = DjangoMigration.objects.filter(app=app)
        max_migration_id = 0

        for migration in migrations:
            # Extract the number from the migration name
            match = re.match(r"^([0-9]+)_", migration.name)  # type: ignore
            if match:
                migration_id = int(match.group(1))
                max_migration_id = max(max_migration_id, migration_id)

        app_migrations[app] = max_migration_id

    current_state = [
        {"migration_id": mig_id, "app": app}
        for app, mig_id in sorted(app_migrations.items())
    ]

    latest_state = get_latest_migration_state()

    if latest_state is None or (current_state != latest_state.migration_state_id):
        MigrationState.objects.create(migration_state_id=current_state)
