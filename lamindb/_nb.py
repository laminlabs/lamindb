from typing import List, Optional, Union

import nbproject as _nb
from lamin_logger import logger
from lnschema_core import Run, Transform

from ._context import context


# this whole class is deprecated, see lamindb.context instead!
class nb:
    """Manage Jupyter notebooks.

    For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.
    """

    transform: Transform = None  # Notebook of this Python session
    run: Run = None  # run of this Python session

    @classmethod
    def header(
        cls,
        *,
        pypackage: Union[str, List[str], None] = None,
        filepath: Optional[str] = None,
        env: Optional[str] = None,
        id: Optional[str] = None,
        v: Optional[str] = "0",
        name: Optional[str] = None,
    ) -> Run:
        """Track the notebook & display metadata.

        Call without arguments in most settings.

        If the notebook has no nbproject metadata, initializes & writes metadata
        to disk.

        Args:
            pypackage: One or more python packages to track.
            filepath: Filepath of notebook. Only needed if automatic inference fails.
            env: Editor environment. Only needed if automatic inference fails.
                Pass `'lab'` for jupyter lab and `'notebook'` for jupyter notebook,
                this can help to identify the correct mechanism for interactivity
                when automatic inference fails.
            id: Pass a notebook id manually.
            v: Pass a notebook version manually.
            name: Pass a notebook name manually.
        """
        logger.warning(
            "DeprecationWarning: Please replace ln.nb.header() with ln.track()."
        )
        context._track(pypackage=pypackage, notebook_path=filepath, editor=env)
        cls.transform = context.transform
        cls.run = context.run
        return cls.run

    @classmethod
    def publish(cls, version: str = None, i_confirm_i_saved: bool = False):
        """Publish the notebook.

        Args:
            version: Notebook version to publish.
            i_confirm_i_saved: Only relevant outside Jupyter Lab & Notebook as a
                safeguard against losing the editor buffer content.
        """
        import lamindb as ln

        if cls.transform is None:
            cls.transform = ln.context.transform
        if cls.run is None:
            cls.run = ln.context.run
        from nbproject._publish import finalize_publish, run_checks_for_publish

        result = run_checks_for_publish(
            calling_statement="publish(", i_confirm_i_saved=i_confirm_i_saved
        )
        if result != "checks-passed":
            return result
        finalize_publish(calling_statement="publish(", version=version)
        # update DB
        transform = ln.select(
            Transform, id=cls.transform.id, version=cls.transform.version  # type: ignore  # noqa
        ).one()
        # update version
        transform.title = _nb.meta.live.title
        if version != cls.transform.version:  # type: ignore
            transform_add = Transform(
                id=transform.id,
                version=version,
                name=transform.name,
                title=transform.title,
            )
        else:
            transform_add = transform
        ln.add(transform_add)
        if version != cls.transform.version:  # type: ignore
            cls.run.transform_version = version  # type: ignore
            ln.add(cls.run)
            ln.delete(transform)
