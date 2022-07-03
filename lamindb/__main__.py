import argparse
from pathlib import Path

from .setup import _setup
from .setup._settings import description

parser = argparse.ArgumentParser(description="Set up lamindb.")
aa = parser.add_argument
aa("command", type=str, choices=["setup"], help="basic setup")
STORAGE_HELP = description.storage_dir
aa("-s", "--storage", type=str, metavar="s", default=None, help=STORAGE_HELP)
USER_HELP = description.user_name
aa("--user", type=str, metavar="s", default=None, help=USER_HELP)
args = parser.parse_args()
INSTANCE_HELP = description.instance_name
aa("--instance", type=str, metavar="s", default=None, help=INSTANCE_HELP)
args = parser.parse_args()


def main():
    if args.command == "setup":
        if args.storage is None:
            storage = input(f"Please paste {description.storage_dir}: ")
        else:
            storage = args.storage
        if args.instance is None:
            if args.storage.startswith(("s3://", "gs://")):
                instance = args.storage.replace("s3://", "")
            else:
                instance = Path(args.storage).stem
        else:
            instance = args.instance
        if args.user is None:
            user = input(f"Please provide your {description.user_name}: ")
        else:
            user = args.user

        _setup.setup_from_cli(
            storage=storage,
            instance=instance,
            user=user,
        )
