from typing import List, Optional, Union

import nbproject as _nb

# from lamin_logger import logger
from lnschema_core import Run, Transform

from ._context import context


# this whole class is deprecated, see lamindb.context instead!
class nb:
    """Manage Jupyter notebooks.

    - Guide: :doc:`guide/run`
    - FAQ: :doc:`/faq/nb`

    For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.
    """

    transform: Transform = None  # Notebook of this Python session
    run: Run = None  # run of this Python session

    @classmethod
    def header(
        cls,
        *,
        run: Optional[str] = None,
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
            run: If `None`, loads the latest run of the notebook (or creates on if
                there isn't any).

                If `"new"`, instantiates a new `Run` record.

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
        # logger.warning(
        #     "DeprecationWarning: Please replace ln.nb.header() with ln.track()."
        # )
        context._track_notebook(
            pypackage=pypackage, filepath=filepath, id=id, v=v, name=name, editor=env
        )
        cls.transform = context.transform
        if run == "new":
            run = Run()
        elif run is None:
            run = Run(load_latest=True)
        else:
            raise ValueError("Pass 'new' to ln.nb.header().")
        cls.run = run
        return run

    @classmethod
    def publish(cls, version: str = None, i_confirm_i_saved: bool = False):
        """Publish the notebook.

        Args:
            version: Notebook version to publish.
            i_confirm_i_saved: Only relevant outside Jupyter Lab & Notebook as a
                safeguard against losing the editor buffer content.
        """
        from nbproject._publish import finalize_publish, run_checks_for_publish

        import lamindb as ln

        result = run_checks_for_publish(
            calling_statement="publish(", i_confirm_i_saved=i_confirm_i_saved
        )
        if result != "checks-passed":
            return result
        finalize_publish(calling_statement="publish(", version=version)
        # update DB
        transform = ln.select(Transform, id=cls.transform.id, v=cls.transform.v).one()
        # update version
        transform.title = _nb.meta.live.title
        if version != cls.transform.v:
            transform_add = Transform(
                id=transform.id, v=version, name=transform.name, title=transform.title
            )
        else:
            transform_add = transform
        ln.add(transform_add)
        if version != cls.transform.v:
            cls.run.transform_v = version
            ln.add(cls.run)
            ln.delete(transform)
