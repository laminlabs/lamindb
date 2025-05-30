# Generated by Django 5.2 on 2025-05-22 15:21

import django.db.models.deletion
import django.db.models.functions.datetime
from django.db import migrations, models

import lamindb.base.fields
import lamindb.base.uids
import lamindb.base.users
import lamindb.models.can_curate
import lamindb.models.run
import lamindb.models.sqlrecord


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0097_remove_schemaparam_param_remove_paramvalue_param_and_more"),
    ]

    operations = [
        migrations.AlterField(
            model_name="feature",
            name="type",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="instances",
                to="lamindb.feature",
            ),
        ),
        migrations.AlterField(
            model_name="project",
            name="type",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="instances",
                to="lamindb.project",
            ),
        ),
        migrations.AlterField(
            model_name="reference",
            name="type",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="instances",
                to="lamindb.reference",
            ),
        ),
        migrations.AlterField(
            model_name="schema",
            name="type",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="instances",
                to="lamindb.schema",
            ),
        ),
        migrations.AlterField(
            model_name="ulabel",
            name="type",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="instances",
                to="lamindb.ulabel",
            ),
        ),
        migrations.CreateModel(
            name="Record",
            fields=[
                (
                    "_branch_code",
                    models.SmallIntegerField(db_default=1, db_index=True, default=1),
                ),
                (
                    "_aux",
                    lamindb.base.fields.JSONField(
                        blank=True, db_default=None, default=None, null=True
                    ),
                ),
                (
                    "created_at",
                    lamindb.base.fields.DateTimeField(
                        blank=True,
                        db_default=django.db.models.functions.datetime.Now(),
                        db_index=True,
                        editable=False,
                    ),
                ),
                (
                    "updated_at",
                    lamindb.base.fields.DateTimeField(
                        blank=True,
                        db_default=django.db.models.functions.datetime.Now(),
                        db_index=True,
                        editable=False,
                    ),
                ),
                ("id", models.AutoField(primary_key=True, serialize=False)),
                (
                    "uid",
                    lamindb.base.fields.CharField(
                        blank=True,
                        db_index=True,
                        default=lamindb.base.uids.base62_16,
                        editable=False,
                        max_length=16,
                        unique=True,
                    ),
                ),
                (
                    "name",
                    lamindb.base.fields.CharField(
                        blank=True,
                        db_index=True,
                        default=None,
                        max_length=150,
                        null=True,
                    ),
                ),
                (
                    "is_type",
                    lamindb.base.fields.BooleanField(
                        blank=True, db_index=True, default=False, null=True
                    ),
                ),
                (
                    "description",
                    lamindb.base.fields.CharField(
                        blank=True, default=None, max_length=255, null=True
                    ),
                ),
                (
                    "created_by",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.base.users.current_user_id,
                        editable=False,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.user",
                    ),
                ),
                (
                    "run",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.models.run.current_run,
                        null=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.run",
                    ),
                ),
                (
                    "space",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        db_default=1,
                        default=1,
                        on_delete=django.db.models.deletion.PROTECT,
                        to="lamindb.space",
                    ),
                ),
                (
                    "type",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="instances",
                        to="lamindb.record",
                    ),
                ),
            ],
            options={
                "abstract": False,
            },
            bases=(lamindb.models.can_curate.CanCurate, models.Model),
        ),
        migrations.CreateModel(
            name="RecordArtifact",
            fields=[
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordartifact",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_artifact",
                        to="lamindb.record",
                    ),
                ),
                (
                    "value",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_record",
                        to="lamindb.artifact",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature", "value")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
        migrations.AddField(
            model_name="record",
            name="artifacts",
            field=models.ManyToManyField(
                related_name="records",
                through="lamindb.RecordArtifact",
                to="lamindb.artifact",
            ),
        ),
        migrations.CreateModel(
            name="RecordProject",
            fields=[
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordproject",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_project",
                        to="lamindb.record",
                    ),
                ),
                (
                    "value",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_record",
                        to="lamindb.project",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
        migrations.AddField(
            model_name="project",
            name="records",
            field=models.ManyToManyField(
                related_name="projects",
                through="lamindb.RecordProject",
                to="lamindb.record",
            ),
        ),
        migrations.CreateModel(
            name="RecordRecord",
            fields=[
                (
                    "_branch_code",
                    models.SmallIntegerField(db_default=1, db_index=True, default=1),
                ),
                (
                    "_aux",
                    lamindb.base.fields.JSONField(
                        blank=True, db_default=None, default=None, null=True
                    ),
                ),
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordrecord",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_record",
                        to="lamindb.record",
                    ),
                ),
                (
                    "space",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        db_default=1,
                        default=1,
                        on_delete=django.db.models.deletion.PROTECT,
                        to="lamindb.space",
                    ),
                ),
                (
                    "value",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_record",
                        to="lamindb.record",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
        migrations.AddField(
            model_name="record",
            name="components",
            field=models.ManyToManyField(
                related_name="composites",
                through="lamindb.RecordRecord",
                to="lamindb.record",
            ),
        ),
        migrations.CreateModel(
            name="RecordRun",
            fields=[
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordrun",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_run",
                        to="lamindb.record",
                    ),
                ),
                (
                    "value",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_record",
                        to="lamindb.run",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
        migrations.AddField(
            model_name="record",
            name="runs",
            field=models.ManyToManyField(
                related_name="records", through="lamindb.RecordRun", to="lamindb.run"
            ),
        ),
        migrations.CreateModel(
            name="RecordULabel",
            fields=[
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordulabel",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_ulabel",
                        to="lamindb.record",
                    ),
                ),
                (
                    "value",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_record",
                        to="lamindb.ulabel",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
        migrations.AddField(
            model_name="record",
            name="ulabels",
            field=models.ManyToManyField(
                related_name="_records",
                through="lamindb.RecordULabel",
                to="lamindb.ulabel",
            ),
        ),
        migrations.CreateModel(
            name="Sheet",
            fields=[
                (
                    "_branch_code",
                    models.SmallIntegerField(db_default=1, db_index=True, default=1),
                ),
                (
                    "_aux",
                    lamindb.base.fields.JSONField(
                        blank=True, db_default=None, default=None, null=True
                    ),
                ),
                (
                    "created_at",
                    lamindb.base.fields.DateTimeField(
                        blank=True,
                        db_default=django.db.models.functions.datetime.Now(),
                        db_index=True,
                        editable=False,
                    ),
                ),
                (
                    "updated_at",
                    lamindb.base.fields.DateTimeField(
                        blank=True,
                        db_default=django.db.models.functions.datetime.Now(),
                        db_index=True,
                        editable=False,
                    ),
                ),
                ("id", models.AutoField(primary_key=True, serialize=False)),
                (
                    "uid",
                    lamindb.base.fields.CharField(
                        blank=True,
                        db_index=True,
                        default=lamindb.base.uids.base62_12,
                        editable=False,
                        max_length=12,
                        unique=True,
                    ),
                ),
                (
                    "name",
                    lamindb.base.fields.CharField(
                        blank=True, db_index=True, default=None, max_length=255
                    ),
                ),
                (
                    "description",
                    lamindb.base.fields.CharField(
                        blank=True,
                        db_index=True,
                        default=None,
                        max_length=255,
                        null=True,
                    ),
                ),
                (
                    "created_by",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.base.users.current_user_id,
                        editable=False,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.user",
                    ),
                ),
                (
                    "run",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.models.run.current_run,
                        null=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.run",
                    ),
                ),
                (
                    "schema",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="sheets",
                        to="lamindb.schema",
                    ),
                ),
                (
                    "space",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        db_default=1,
                        default=1,
                        on_delete=django.db.models.deletion.PROTECT,
                        to="lamindb.space",
                    ),
                ),
            ],
            options={
                "abstract": False,
            },
        ),
        migrations.AddField(
            model_name="record",
            name="sheet",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="records",
                to="lamindb.sheet",
            ),
        ),
        migrations.CreateModel(
            name="SheetProject",
            fields=[
                (
                    "created_at",
                    lamindb.base.fields.DateTimeField(
                        blank=True,
                        db_default=django.db.models.functions.datetime.Now(),
                        db_index=True,
                        editable=False,
                    ),
                ),
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "created_by",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.base.users.current_user_id,
                        editable=False,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.user",
                    ),
                ),
                (
                    "project",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="links_sheet",
                        to="lamindb.project",
                    ),
                ),
                (
                    "run",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        default=lamindb.models.run.current_run,
                        null=True,
                        on_delete=django.db.models.deletion.PROTECT,
                        related_name="+",
                        to="lamindb.run",
                    ),
                ),
                (
                    "sheet",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_project",
                        to="lamindb.sheet",
                    ),
                ),
            ],
            options={
                "unique_together": {("sheet", "project")},
            },
            bases=(lamindb.models.sqlrecord.IsLink, models.Model),
        ),
        migrations.AddField(
            model_name="project",
            name="sheets",
            field=models.ManyToManyField(
                related_name="projects",
                through="lamindb.SheetProject",
                to="lamindb.sheet",
            ),
        ),
        migrations.CreateModel(
            name="RecordJson",
            fields=[
                ("id", models.BigAutoField(primary_key=True, serialize=False)),
                (
                    "value",
                    lamindb.base.fields.JSONField(
                        blank=True, db_default=None, default=None
                    ),
                ),
                (
                    "feature",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="links_recordjson",
                        to="lamindb.feature",
                    ),
                ),
                (
                    "record",
                    lamindb.base.fields.ForeignKey(
                        blank=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="values_json",
                        to="lamindb.record",
                    ),
                ),
            ],
            options={
                "unique_together": {("record", "feature")},
            },
            bases=(models.Model, lamindb.models.sqlrecord.IsLink),
        ),
    ]
