from pathlib import Path

import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import dobject

from lamindb import schema


def storage_key_from_dobject(dobj: dobject):
    return f"{dobj.id}-{dobj.v}{dobj.file_suffix}"


def storage_key_from_triple(dobj_id: str, dobj_v: str, dobj_suffix: str):
    return f"{dobj_id}-{dobj_v}{dobj_suffix}"


def filepath_from_dobject(dobj: dobject):
    storage_key = storage_key_from_dobject(dobj)
    filepath = settings.instance.storage.key_to_filepath(storage_key)
    return filepath


def track_usage(dobject_id, dobject_v, usage_type):
    with sqm.Session(settings.instance.db_engine()) as session:
        usage = schema.core.usage(
            type=usage_type,
            user_id=settings.user.id,
            dobject_id=dobject_id,
            dobject_v=dobject_v,
        )
        session.add(usage)
        session.commit()
        session.refresh(usage)

    settings.instance._update_cloud_sqlite_file()

    return usage.id


def format_pipeline_logs(logs):
    pipeline_dir_logs = {}
    for log in logs:
        rel_filepath = Path(log[0].split(" ")[0])
        print(rel_filepath)
        n_parents = len(rel_filepath.parents) - 1
        if n_parents != 0:
            logs.remove(log)
            top_dir = str(list(rel_filepath.parents)[-2])
            if top_dir in pipeline_dir_logs.keys():
                pipeline_dir_logs[top_dir][0] += 1
            else:
                jupynb_log = log[1]
                settings_log = log[2]
                pipeline_dir_logs[top_dir] = (1, jupynb_log, settings_log)

    for dir, logs in pipeline_dir_logs.items():
        logs.append([f"{logs[0]} files in {dir}", logs[1], logs[2]])

    return logs
