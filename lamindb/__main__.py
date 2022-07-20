import argparse
from pathlib import Path

from ._setup._settings import description
from .refactoring.Auth import Auth
from .refactoring.SetupInstance import SetupInstance

description_cli = """
Configure LaminDB and perform simple actions with these commands:
- lndb signup --email <email> | First time sign up & log in after email is confirmed.
- lndb login [--email <email>] [--secret <secret>] | Log in an already-signed-up user.
- lndb logout | Log in an already-signed-up user.
- lndb init [--name <name>] | Init instance.
- lndb remove [--name <name>] | Remove instance.
"""  # noqa
parser = argparse.ArgumentParser(
    description=description_cli, formatter_class=argparse.RawTextHelpFormatter
)
aa = parser.add_argument
aa("command", type=str, choices=["signup", "login", "logout", "init", "remove"])
# user
aa("--email", type=str, metavar="s", default=None, help=description.user_email)
aa("--secret", type=str, metavar="s", default=None, help=description.user_secret)
# db instance
aa("--name", type=str, metavar="name", default=None, help=description.instance_name)
aa("--storage", type=str, metavar="s", default=None, help=description.storage_dir)
aa("--db", type=str, metavar="s", default="sqlite", help=description._dbconfig)
args = parser.parse_args()


def main():
    auth = Auth()

    if args.command == "signup":
        return auth.sign_up(
            email=args.email,
            secret=args.secret,
        )
    elif args.command == "login":
        return auth.log_in(
            email=args.email,
            secret=args.secret,
        )
    elif args.command == "logout":
        return auth.log_out()
    elif args.command == "init":
        return SetupInstance.setup_if_not_exists(
            instance_name=args.name,
            db_base_path=Path(f"./{args.name}/db"),
            storage_base_path=Path(f"./{args.name}/storage"),
        )
    elif args.command == "remove":
        return SetupInstance.remove_instance(instance_name=args.name)
    else:
        raise RuntimeError("Invalid command. Allowed are: `lndb user` & `lndb db`")
