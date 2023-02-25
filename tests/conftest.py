from pathlib import Path
from subprocess import run


def pytest_sessionstart(session):
    instance_dirs = [
        d for d in ["./docs/guide/mydata", "./mydata-test-db"] if Path(d).exists()
    ]
    for instance_dir in instance_dirs:
        clean_instance = f"rm -r {instance_dir}"
        session.run(*clean_instance.split(" "))
    run("lamin init --storage mydata-test-db", shell=True)
