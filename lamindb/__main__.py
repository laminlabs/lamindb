import argparse

from .setup import setup
from .setup._settings import description

parser = argparse.ArgumentParser(description="Set up lamindb.")
aa = parser.add_argument
aa("command", type=str, choices=["setup"], help="basic setup")
STORAGE_HELP = description.storage_root
aa("-s", "--storage", type=str, metavar="s", default=None, help=STORAGE_HELP)
CACHE_HELP = description.cache_root
aa("--cache", type=str, metavar="s", default=None, help=CACHE_HELP)
USER_HELP = description.user_name
aa("--user", type=str, metavar="s", default=None, help=USER_HELP)
args = parser.parse_args()


def main():
    if args.command == "setup":
        setup(
            storage=args.storage,
            cache=args.cache,
            user=args.user,
        )
