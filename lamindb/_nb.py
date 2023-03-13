from typing import List, Optional, Union

import nbproject as _nb
from lnschema_core import Notebook, Run

from ._context import context


# this whole class is deprecated, see lamindb.context instead!
class nb:
    """Manage Jupyter notebooks.

    - Guide: :doc:`guide/run`
    - FAQ: :doc:`/faq/nb`

    For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.
    """

    notebook: Notebook = None  # Notebook of this Python session
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
    ):
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
        context._track_notebook(
            pypackage=pypackage, filepath=filepath, id=id, v=v, name=name, editor=env
        )
        notebook = context.notebook
        cls.notebook = notebook
        if run == "new":
            Run(global_context=True)
        elif run is None:
            Run(global_context=True, load_latest=True)
        else:
            raise ValueError("Pass 'new' to ln.nb.header().")

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
        import lamindb.schema as lns

        result = run_checks_for_publish(
            calling_statement="publish(", i_confirm_i_saved=i_confirm_i_saved
        )
        if result != "checks-passed":
            return result
        finalize_publish(calling_statement="publish(", version=version)
        # update DB
        notebook = ln.select(Notebook, id=cls.notebook.id, v=cls.notebook.v).one()
        # update version
        notebook.title = _nb.meta.live.title
        if version != cls.notebook.v:
            notebook_add = lns.Notebook(
                id=notebook.id, v=version, name=notebook.name, title=notebook.title
            )
        else:
            notebook_add = notebook
        ln.add(notebook_add)
        if version != cls.notebook.v:
            cls.run.notebook_v = version
            ln.add(cls.run)
            ln.delete(notebook)
