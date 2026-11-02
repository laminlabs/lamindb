import django.db.models.deletion
from django.db import migrations, models

import lamindb.models.run


class Migration(migrations.Migration):
    dependencies = [("lamindb", "0193_v2_9_part_2")]

    operations = [
        migrations.AddField(
            model_name="recordblock",
            name="previous_version",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.RESTRICT,
                related_name="next_versions",
                to="lamindb.recordblock",
            ),
        ),
        migrations.AddField(
            model_name="recordblock",
            name="run",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="output_record_blocks",
                to="lamindb.run",
            ),
        ),
        # Add the Python default only after the nullable column exists: migration
        # must not attribute historical rows to any active run.
        migrations.AlterField(
            model_name="recordblock",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.run.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="output_record_blocks",
                to="lamindb.run",
            ),
        ),
    ]
