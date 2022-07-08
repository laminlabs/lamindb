import argparse

from ._setup import _setup
from ._setup._settings import description

description_cli = """
Configure LaminDB and perform simple actions with these commands:
- lndb signup --email <email> | First time sign up & log in after email is confirmed.
- lndb login [--email <email>] [--secret <secret>] | Log in an already-signed-up user.
- lndb config [--storage <storage>] [--db <db>] | Configure/switch storage or instance.
"""
parser = argparse.ArgumentParser(
    description=description_cli, formatter_class=argparse.RawTextHelpFormatter
)
aa = parser.add_argument
aa("command", type=str, choices=["signup", "login", "config"])
# user
aa("--email", type=str, metavar="s", default=None, help=description.user_email)
aa("--secret", type=str, metavar="s", default=None, help=description.user_secret)
# db instance
aa("--storage", type=str, metavar="s", default=None, help=description.storage_dir)
aa("--db", type=str, metavar="s", default="sqlite", help=description._dbconfig)
args = parser.parse_args()


def main():
    if args.command == "signup":
        return _setup.sign_up_first_time(
            email=args.email,
        )
    if args.command == "login":
        return _setup.log_in_user(
            email=args.email,
            secret=args.secret,
        )
    elif args.command == "config":
        return _setup.setup_instance(
            storage=args.storage,
            dbconfig=args.db,
        )
    else:
        raise RuntimeError("Invalid command. Allowed are: `lndb user` & `lndb db`")
