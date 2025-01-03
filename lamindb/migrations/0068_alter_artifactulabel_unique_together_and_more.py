# Generated by Django 5.1.1 on 2024-10-18 14:31

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0067_alter_featurevalue_unique_together_and_more"),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name="artifactulabel",
            unique_together={("artifact", "ulabel", "feature")},
        ),
        migrations.AlterUniqueTogether(
            name="featuresetfeature",
            unique_together={("featureset", "feature")},
        ),
    ]
