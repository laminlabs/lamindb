import django.db.models.functions.text
from django.db import migrations, models


def lowercase_default_values(apps, schema_editor):
    """Lowercase default values for Space and Branch models."""
    Space = apps.get_model("lamindb", "Space")
    Branch = apps.get_model("lamindb", "Branch")

    space = Space.objects.get(uid="A")
    space.uid = "a"
    space.name = "all"
    space.save()

    trash_branch = Branch.objects.get(uid="T")
    trash_branch.uid = "t"
    trash_branch.name = "trash"
    trash_branch.save()

    archive_branch = Branch.objects.get(uid="A")
    archive_branch.uid = "a"
    archive_branch.name = "archive"
    archive_branch.save()

    main_branch = Branch.objects.get(uid="M")
    main_branch.uid = "m"
    main_branch.name = "main"
    main_branch.save()


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0112_alter_recordartifact_feature_and_more"),
    ]

    operations = [
        migrations.RunPython(
            lowercase_default_values,
        ),
        migrations.AlterModelOptions(
            name="branch",
            options={},
        ),
        migrations.AlterModelOptions(
            name="space",
            options={},
        ),
        migrations.AddConstraint(
            model_name="branch",
            constraint=models.UniqueConstraint(
                django.db.models.functions.text.Lower("name"),
                name="unique_branch_name_lower",
            ),
        ),
        migrations.AddConstraint(
            model_name="space",
            constraint=models.UniqueConstraint(
                django.db.models.functions.text.Lower("name"),
                name="unique_space_name_lower",
            ),
        ),
    ]
