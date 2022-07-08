import argparse

from ._setup import _setup
from ._setup._settings import description

parser = argparse.ArgumentParser(description="Set up lamindb.")
aa = parser.add_argument
aa("command", type=str, choices=["setup"], help="basic setup")
STORAGE_HELP = description.storage_dir
aa("-s", "--storage", type=str, metavar="s", default=None, help=STORAGE_HELP)
USER_HELP = description.user_email
aa("--user", type=str, metavar="s", default=None, help=USER_HELP)
SECRET_HELP = description.user_secret
aa("--secret", type=str, metavar="s", default=None, help=SECRET_HELP)
INSTANCE_HELP = description.instance_name
aa("--instance", type=str, metavar="s", default=None, help=INSTANCE_HELP)
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
        )
