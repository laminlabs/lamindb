# Generated by Django 5.2 on 2024-08-02 08:03

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0057_link_models_latest_report_and_others"),
    ]

    operations = [
        migrations.AddField(
            model_name="artifact",
            name="_actions",
            field=models.ManyToManyField(related_name="+", to="lamindb.artifact"),
        ),
        migrations.AddField(
            model_name="collection",
            name="_actions",
            field=models.ManyToManyField(related_name="+", to="lamindb.artifact"),
        ),
    ]
