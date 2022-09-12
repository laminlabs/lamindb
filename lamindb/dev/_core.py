from pathlib import Path

import sqlmodel as sqm
from lndb_setup import settings
from lnschema_core import dobject, type

from lamindb import schema


def storage_key_from_dobject(dobj: dobject):
    return f"{dobj.id}-{dobj.v}{dobj.file_suffix}"


def storage_key_from_triple(dobj_id: str, dobj_v: str, dobj_suffix: str):
    return f"{dobj_id}-{dobj_v}{dobj_suffix}"


def filepath_from_dobject(dobj: dobject):
    storage_key = storage_key_from_dobject(dobj)
    filepath = settings.instance.storage.key_to_filepath(storage_key)
    return filepath


def track_usage(dobject_id, dobject_v, usage_type: type.usage):
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
    from operator import itemgetter

    run_logs = {}
    for log in list(logs):
        # parse number of parent directories within run directory
        rel_filepath = Path(log[0].split(" ")[0])
        n_parents = len(rel_filepath.parents) - 1
        if n_parents > 1:
            # remove log and register in dictionnary under its top directory
            logs.remove(log)
            run_log = log[1]
            run_id = run_log.split(" ")[1]
            top_dir = str(list(rel_filepath.parents)[-3])
            if run_id not in run_logs.keys():
                run_logs[run_id] = {}
            dir_logs = run_logs[run_id]
            if top_dir in dir_logs.keys():
                dir_logs[top_dir][0] += 1
            else:
                run_dir_log = log[2]
                user_log = log[3]
                dir_logs[top_dir] = [1, run_log, run_dir_log, user_log]

    run_dir_logs = [
        (dir, logs) for run in run_logs.values() for (dir, logs) in run.items()
    ]

    for dir, new_logs in run_dir_logs:
        logs.append(
            [f"{new_logs[0]} files in {dir}/", new_logs[1], new_logs[2], new_logs[3]]
        )

    logs = sorted(logs, key=itemgetter(1))

    return logs
