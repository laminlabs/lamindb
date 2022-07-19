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
- lndb init [--storage <storage>] [--db <db>] | Init & config instance or storage.
"""  # noqa
parser = argparse.ArgumentParser(
    description=description_cli, formatter_class=argparse.RawTextHelpFormatter
)
aa = parser.add_argument
aa("command", type=str, choices=["signup", "login", "logout", "init"])
# user
aa("--email", type=str, metavar="s", default=None, help=description.user_email)
aa("--secret", type=str, metavar="s", default=None, help=description.user_secret)
# db instance
aa("--storage", type=str, metavar="s", default=None, help=description.storage_dir)
aa("--db", type=str, metavar="s", default="sqlite", help=description._dbconfig)
args = parser.parse_args()


def main():
    supabase_client_url = "https://qntuvxhregqtypdorhaw.supabase.co"
    supabase_client_anon_key = (
        "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9"
        ".eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InFudHV2eGhyZWdxdHlwZG9yaGF3Iiw"
        "icm9sZSI6ImFub24iLCJpYXQiOjE2NTcyMjMxODksImV4cCI6MTk3Mjc5OTE4OX0"
        ".hzrRK-xbhFSqFb1-cVu6NFM-UcPON4HBdaH2qe3vKbA"
    )

    auth = Auth(supabase_client_url, supabase_client_anon_key)

    instance_name = "abc"
    settings_base_path = Path(
        "/Users/fredericenard/Sources/lamin/lamindb/lamindb_data/settings"
    )
    db_base_path = Path("/Users/fredericenard/Sources/lamin/lamindb/lamindb_data/db")
    storage_base_path = Path(
        "/Users/fredericenard/Sources/lamin/lamindb/lamindb_data/storage"
    )

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
            instance_name=instance_name,
            settings_base_path=settings_base_path,
            db_base_path=db_base_path,
            storage_base_path=storage_base_path,
        )
    else:
        raise RuntimeError("Invalid command. Allowed are: `lndb user` & `lndb db`")
