from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

from ..dev._docs import doc_args
from ._settings import Settings, description
from ._setup_db import setup_db


def setup_storage(
    storage_root: Union[str, Path, CloudPath] = None,
    cache_root: Union[str, Path] = None,
):

    if storage_root is None:
        storage_root = input(f"Please paste {description.storage_root}: ")

    # check whether a local directory actually exists
    if isinstance(storage_root, str) and storage_root.startswith(("s3://", "gs://")):
        cloud_storage = True
        storage_root = storage_root
    else:
        cloud_storage = False
        storage_root = Path(storage_root)
        if not storage_root.exists():
            print(f"creating storage directory {storage_root}")
            storage_root.mkdir(parents=True)

    if cloud_storage:
        # define cache directory
        parents = [Path("/users/shared/"), Path.home()]
        examples = ", or ".join([str(p / ".lamindb" / "cache") for p in parents])
        if cache_root is None:
            cache_root = input(
                f"Please paste {description.cache_root}, e.g., {examples}: "
            )
        cache_root = Path(cache_root)
        if not cache_root.exists():
            print(f"creating cache directory {cache_root}")
            cache_root.mkdir(parents=True)
        else:
            print(f"using cache directory {cache_root}")
    else:
        # we do not need a cache as we're not working in the cloud
        cache_root = None

    return dict(
        storage_root_str=str(storage_root),
        cache_root=cache_root,
    )


def setup_user(user: str = None):

    if user is None:
        user = input(f"Please provide your {description.user_name}: ")

    return dict(user_name=user)


def setup_notion(notion: str = None):

    NOTION_HELP = (
        "Notion integration token (currently dysfunctional, see"
        " https://lamin.ai/notes/2022/explore-notion)"
    )

    if notion is None:
        notion = input(
            f"Please paste your internal {NOTION_HELP}"
            " (https://notion.so/my-integrations): "
        )

    return dict(NOTION_API_KEY=notion)


@doc_args(description.storage_root, description.cache_root, description.user_name)
def setup(
    *,
    storage: str = None,
    cache: str = None,
    user: str = None,
):
    """Setup LaminDB. Alternative to using the CLI via `lamindb setup`.

    Args:
        storage: {}.
        cache: {}.
        user: {}.
    """
    settings = {}
    settings.update(setup_storage(storage, cache))
    settings.update(setup_user(user))

    Settings(**settings)._write()  # type: ignore

    setup_db(settings)  # update settings with user_id

    settings = Settings(**settings)  # type: ignore

    settings._write()  # type: ignore

    print("Successfully set up lamindb!")

    return settings
