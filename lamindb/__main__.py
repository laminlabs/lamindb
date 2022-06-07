import argparse
from pathlib import Path
from typing import Union

# consider changing to click
# * https://click.palletsprojects.com/en/5.x/
# * https://collectiveacuity.medium.com/argparse-vs-click-227f53f023dc
parser = argparse.ArgumentParser(description="Setup laminlake.")
aa = parser.add_argument
aa("command", type=str, choices=["configure"], help="basic setup")
STORAGE_HELP = (
    "storage root, if not a local directory, it needs to be of form 's3://bucket_name'"
    " or 'gs://bucket_name'"
)
aa("-s", "--storage", type=str, metavar="s", default=None, help=STORAGE_HELP)
CACHE_HELP = "cache root, a local directory to cache cloud files"
aa("--cache", type=str, metavar="s", default=None, help=STORAGE_HELP)
NOTION_HELP = "Notion integration token"
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
        examples = ", or ".join(
            [
                str(p / ".lamin" / "cache")
                for p in [
                    Path("/users/shared/"),
                    Path.home(),
                ]
            ]
        )
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


def configure_notion(notion: str = None):

    if notion is None:
        notion = input(
            f"Please paste your internal {NOTION_HELP}"
            " (https://notion.so/my-integrations): "
        )

    # write a _secrets.py file that's in .gitignore
    with open(root_dir / "_secrets.py", "w") as f:
        f.write(f"NOTION_API_KEY = {notion!r}\n")


def main():
    if args.command == "configure":
        configure_storage(storage_root=args.storage, cache_root=args.cache)
        if args.notion is not None:
            configure_notion(notion=args.notion)
