# Generated by Django 5.2 on 2024-09-09 07:52

import django.db.models.deletion
from django.db import migrations, models

import lamindb.base.users
import lamindb.models


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0065_remove_collection_feature_sets_and_more"),
    ]

    operations = [
        migrations.AlterField(
            model_name="artifact",
            name="_feature_values",
            field=models.ManyToManyField(
                related_name="artifacts",
                through="lamindb.ArtifactFeatureValue",
                to="lamindb.featurevalue",
            ),
        ),
        migrations.AlterField(
            model_name="artifact",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="artifactfeatureset",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="artifactfeatureset",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="artifactfeaturevalue",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="artifactfeaturevalue",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="artifactulabel",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="artifactulabel",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="collection",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="collectionartifact",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="collectionartifact",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="collectionulabel",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="collectionulabel",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="feature",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="feature",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="featureset",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="featureset",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="featurevalue",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="featurevalue",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="param",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="param",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="storage",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="storage",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="ulabel",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="ulabel",
            name="run",
            field=models.ForeignKey(
                default=lamindb.models.current_run,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="artifact",
            name="_previous_runs",
            field=models.ManyToManyField(
                related_name="_output_artifacts_with_later_updates",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="artifact",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="created_artifacts",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="collection",
            name="_previous_runs",
            field=models.ManyToManyField(
                related_name="_output_collections_with_later_updates",
                to="lamindb.run",
            ),
        ),
        migrations.AlterField(
            model_name="featurevalue",
            name="feature",
            field=models.ForeignKey(
                default=None,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="values",
                to="lamindb.feature",
            ),
        ),
        migrations.AlterField(
            model_name="paramvalue",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="+",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="paramvalue",
            name="param",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="values",
                to="lamindb.param",
            ),
        ),
        migrations.AlterField(
            model_name="run",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="created_runs",
                to="lamindb.user",
            ),
        ),
        migrations.AlterField(
            model_name="transform",
            name="created_by",
            field=models.ForeignKey(
                default=lamindb.base.users.current_user_id,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="created_transforms",
                to="lamindb.user",
            ),
        ),
    ]
