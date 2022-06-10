import argparse
from pathlib import Path
from typing import Union

# consider changing to click
# * https://click.palletsprojects.com/en/5.x/
# * https://collectiveacuity.medium.com/argparse-vs-click-227f53f023dc
parser = argparse.ArgumentParser(description="Set up lamindb.")
aa = parser.add_argument
aa("command", type=str, choices=["setup"], help="basic setup")
STORAGE_HELP = (
    "storage root, if not a local directory, it needs to be of form 's3://bucket_name'"
    " or 'gs://bucket_name'"
)
aa("-s", "--storage", type=str, metavar="s", default=None, help=STORAGE_HELP)
CACHE_HELP = "cache root, a local directory to cache cloud files"
aa("--cache", type=str, metavar="s", default=None, help=CACHE_HELP)
USER_HELP = "(GitHub) user name"
aa("--user", type=str, metavar="s", default=None, help=USER_HELP)
NOTION_HELP = (
    "Notion integration token (currently dysfunctional, see"
    " https://lamin.ai/notes/2022/explore-notion)"
)
aa("--notion", type=str, metavar="token", default=None, help=NOTION_HELP)
args = parser.parse_args()

root_dir = Path(__file__).parent.resolve()


def configure_storage(
    storage_root: Union[str, Path] = None, cache_root: Union[str, Path] = None
):

    if storage_root is None:
        storage_root = input(f"Please paste {STORAGE_HELP}: ")

    # check whether a local directory actually exists
    if isinstance(storage_root, str) and storage_root.startswith(("s3://", "gs://")):
        cloud_storage = True
    else:
        cloud_storage = False
        storage_root = Path(storage_root)
        if not storage_root.exists():
            print(f"creating storage directory {storage_root}")
            storage_root.mkdir(parents=True)

    if cloud_storage:
        # define cache directory
        parents = [Path("/users/shared/"), Path.home()]
        examples = ", or ".join([str(p / ".lamin" / "cache") for p in parents])
        if cache_root is None:
            cache_root = input(f"Please paste {CACHE_HELP}, e.g., {examples}: ")
        cache_root = Path(cache_root)
        if not cache_root.exists():
            print(f"creating cache directory {cache_root}")
            cache_root.mkdir(parents=True)
    else:
        # we do not need a cache as we're not working in the cloud
        cache_root = ""

    with open(root_dir / "_configuration.py", "w") as f:
        f.write(f"cloud_storage = {cloud_storage}\n")
        f.write(f"storage_root = {str(storage_root)!r}\n")
        f.write(f"cache_root = {str(cache_root)!r}\n")


def configure_user(user: str = None):

    if user is None:
        user = input(f"Please provide your {USER_HELP}: ")

    # write a _secrets.py file that's in .gitignore
    with open(root_dir / "_configuration.py", "a") as f:
        f.write(f"user_name = {user!r}\n")


def configure_notion(notion: str = None):

    if notion is None:
        notion = input(
            f"Please paste your internal {NOTION_HELP}"
            " (https://notion.so/my-integrations): "
        )

    # write a _secrets.py file that's in .gitignore
    with open(root_dir / "_secrets.py", "w") as f:
        f.write(f"NOTION_API_KEY = {notion!r}\n")


def setup():
    configure_storage(storage_root=args.storage, cache_root=args.cache)
    configure_user(user=args.user)
    if args.notion is not None:
        configure_notion(notion=args.notion)

    # set up database
    from lamindb import db

    db.meta.create()

    from lamindb._configuration import user_name

    user_id = db.insert_if_not_exists.user(user_name)

    # write a _secrets.py file that's in .gitignore
    with open(root_dir / "_configuration.py", "a") as f:
        f.write(f"user_id = {user_id!r}\n")

    print("successfully set up lamindb!")


def main():
    if args.command == "setup":
        setup()
