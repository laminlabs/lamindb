import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [("lamindb", "0193_v2_9_part_2")]

    operations = [
        migrations.CreateModel(
            name="ArtifactUIDReservation",
            fields=[
                (
                    "uid",
                    models.CharField(max_length=20, primary_key=True, serialize=False),
                ),
                ("key", models.CharField(max_length=1024, null=True)),
                ("created_at", models.DateTimeField(auto_now_add=True)),
                ("state", models.CharField(default="reserved", max_length=20)),
                (
                    "storage",
                    models.ForeignKey(
                        null=True,
                        on_delete=django.db.models.deletion.SET_NULL,
                        to="lamindb.storage",
                    ),
                ),
            ],
        ),
    ]
