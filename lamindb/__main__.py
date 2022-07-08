import argparse

from ._setup import _setup
from ._setup._settings import description

parser = argparse.ArgumentParser(description="Configure LaminDB.")
aa = parser.add_argument
aa("command", type=str, choices=["set"], help="Set settings.")
aa("--user", type=str, metavar="s", default=None, help=description.user_email)
aa("--secret", type=str, metavar="s", default=None, help=description.user_secret)
aa("--storage", type=str, metavar="s", default=None, help=description.storage_dir)
aa("--instance", type=str, metavar="s", default=None, help=description.instance_name)
aa("--db", type=str, choices=["set"], default="sqlite", help=description.db)
args = parser.parse_args()


def main():
    if args.command == "setup":
        if args.storage is None:
            storage = input(f"Please paste {description.storage_dir}: ")
        else:
            storage = args.storage

        return _setup.setup_from_cli(
            storage=storage,
            user=args.user,
            secret=args.secret,
            instance=args.instance,
            db=args.instance,
        )
