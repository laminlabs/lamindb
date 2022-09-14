from pathlib import Path


def format_pipeline_logs(logs):
    from operator import itemgetter

    run_logs = {}
    for log in list(logs):
        # parse number of parent directories within run directory
        rel_filepath = Path(log[0].split(" ")[0])
        n_parents = len(rel_filepath.parents)
        if n_parents > 1:
            # remove log and register in dictionnary under its top directory
            logs.remove(log)
            run_log = log[1]
            run_id = run_log.split(" ")[1]
            top_dir = str(list(rel_filepath.parents)[-2])
            if run_id not in run_logs.keys():
                run_logs[run_id] = {}
            dir_logs = run_logs[run_id]
            if top_dir in dir_logs.keys():
                dir_logs[top_dir][0] += 1
            else:
                user_log = log[2]
                dir_logs[top_dir] = [1, run_log, user_log]

    run_dir_logs = [
        (dir, logs) for run in run_logs.values() for (dir, logs) in run.items()
    ]

    for dir, new_logs in run_dir_logs:
        logs.append([f"{new_logs[0]} files in {dir}/", new_logs[1], new_logs[2]])

    logs = sorted(logs, key=itemgetter(1))

    return logs
