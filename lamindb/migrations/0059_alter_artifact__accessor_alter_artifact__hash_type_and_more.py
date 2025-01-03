# Generated by Django 5.2 on 2024-08-02 15:03

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0058_artifact__actions_collection__actions"),
    ]

    operations = [
        migrations.AlterField(
            model_name="artifact",
            name="_accessor",
            field=models.CharField(
                db_index=True, default=None, max_length=64, null=True
            ),
        ),
        migrations.AlterField(
            model_name="artifact",
            name="_hash_type",
            field=models.CharField(
                db_index=True, default=None, max_length=30, null=True
            ),
        ),
        migrations.AlterField(
            model_name="artifact",
            name="_key_is_virtual",
            field=models.BooleanField(),
        ),
    ]
