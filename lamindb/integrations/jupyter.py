"""Notebook helpers for tracking and finishing Jupyter notebooks.

Adapted from nbproject (same authors). This keeps only what lamindb calls:
locating the running notebook and reading its first markdown title.
Notebooks are read with nbformat.
"""

from __future__ import annotations

import json
import os
import sys
from itertools import chain
from pathlib import Path, PurePath
from urllib import request

from lamin_utils import logger

_DIR_KEYS = ("notebook_dir", "root_dir")
_CONN_ERROR = (
    "Unable to access server;\n"
    "querying requires either no security or token based security."
)


def read_notebook(filepath: str | Path):
    """Read a notebook from disk.

    Raises:
        ImportError: If `nbformat` or `jupytext` is not installed.
    """
    try:
        import jupytext  # noqa: F401
        import nbformat
    except ImportError as error:
        raise ImportError(
            "install nbconvert & jupytext: pip install nbconvert jupytext"
        ) from error

    with Path(filepath).open(encoding="utf-8") as file:
        return nbformat.read(file, as_version=4)


def _cell_source(cell) -> str:
    source = cell["source"]
    if isinstance(source, str):
        return source
    return "".join(source)


def get_title(nb) -> str | None:
    """Return the first level-1 markdown heading."""
    for cell in nb.cells:
        if cell["cell_type"] != "markdown":
            continue
        for line in _cell_source(cell).split("\n"):
            if line.startswith("# "):
                return line.lstrip("#").strip(" .").strip("\n")
    return None


def _prepare_url(server: dict, query_str: str = "") -> str:
    token = server["token"]
    if token:
        query_str = f"{query_str}?token={token}"
    return f"{server['url']}api/sessions{query_str}"


def _query_server(server: dict):
    try:
        url = _prepare_url(server)
        with request.urlopen(url) as req:  # noqa: S310
            return json.loads(req.read())
    except Exception:
        raise Exception(_CONN_ERROR) from None


def _running_servers() -> tuple[list, list]:
    nbapp_import = False
    try:
        from notebook.notebookapp import list_running_servers

        nbapp_import = True
        servers_nbapp = list(list_running_servers())
    except ModuleNotFoundError:
        servers_nbapp = []

    try:
        from jupyter_server.serverapp import list_running_servers

        servers_juserv = list(list_running_servers())
    except ModuleNotFoundError:
        servers_juserv = []
        if not nbapp_import:
            logger.warning(
                "It looks like you are running jupyter lab "
                "but don't have jupyter-server module installed. "
                "Please install it via pip install jupyter-server"
            )
    return servers_nbapp, servers_juserv


def _ipylab_is_installed() -> bool:
    from importlib.util import find_spec

    return find_spec("ipylab") is not None


def _lab_notebook_path() -> Path | None:
    try:
        from ipylab import JupyterFrontEnd
    except ImportError:
        return None
    try:
        current_session = JupyterFrontEnd().sessions.current_session
    except Exception:
        return None
    if not current_session or "name" not in current_session:
        return None
    return Path.cwd() / current_session["name"]


def _find_nb_path_via_parent_process() -> Path | None:
    """Find the notebook path from the parent process command line.

    Used when a notebook is executed with nbconvert. Requires psutil.
    """
    import psutil

    try:
        current_process = psutil.Process(os.getpid())
        parent_process = current_process.parent()
        if parent_process is None:
            logger.warning("psutil: Could not get parent process.")
            return None

        cmdline = parent_process.cmdline()
        if not cmdline:
            logger.warning(
                f"psutil: Parent process ({parent_process.pid}) has empty cmdline."
            )
            return None

        logger.info(f"psutil: Parent cmdline: {cmdline}")

        is_nbconvert_call = False
        potential_path = None
        for i, arg in enumerate(cmdline):
            if "nbconvert" in arg.lower():
                base_arg = Path(arg).name.lower()
                if (
                    "jupyter-nbconvert" in base_arg
                    or arg == "nbconvert"
                    or (
                        cmdline[i - 1].endswith("python")
                        and arg == "-m"
                        and cmdline[i + 1] == "nbconvert"
                    )
                ):
                    is_nbconvert_call = True
            if arg.endswith(".ipynb"):
                potential_path = arg

        if is_nbconvert_call and "--inplace" not in cmdline:
            raise ValueError(
                "Please execute notebook 'nbconvert' by passing option '--inplace'."
            )

        if is_nbconvert_call and potential_path:
            try:
                parent_cwd = parent_process.cwd()
                resolved_path = Path(parent_cwd) / Path(potential_path)
                if resolved_path.is_file():
                    logger.info(f"psutil: Found potential path: {resolved_path}")
                    return resolved_path.resolve()
                abs_path = Path(potential_path)
                if abs_path.is_absolute() and abs_path.is_file():
                    logger.info(f"psutil: Found potential absolute path: {abs_path}")
                    return abs_path.resolve()
                logger.warning(
                    f"psutil: Potential path '{potential_path}' not found relative to"
                    f" parent CWD '{parent_cwd}' or as absolute path."
                )
                return None
            except psutil.AccessDenied:
                logger.warning("psutil: Access denied when getting parent CWD.")
                maybe_path = Path(potential_path)
                if maybe_path.is_file():
                    return maybe_path.resolve()
                return None
            except Exception as error:
                logger.warning(
                    f"psutil: Error resolving path '{potential_path}': {error}"
                )
                return None

        logger.warning(
            "psutil: Could not reliably identify notebook path from parent cmdline."
        )
        return None
    except ImportError:
        logger.warning("psutil library not found. Cannot inspect parent process.")
        return None
    except psutil.Error as error:
        logger.warning(f"psutil error: {error}")
        return None
    except ValueError:
        raise
    except Exception as error:
        logger.warning(f"Unexpected error during psutil check: {error}")
        return None


def _with_env(nb_path, env: str | None, return_env: bool, default_env: str):
    if return_env:
        return nb_path, default_env if env is None else env
    return nb_path


def notebook_path(return_env: bool = False):
    """Return the path to the current notebook.

    Args:
        return_env: If `True`, also return where the notebook is running:
            `'lab'`, `'notebook'`, `'vs_code'`, `'nbconvert'`, or `'test'`.
    """
    env = os.environ.get("NBPRJ_TEST_NBENV")
    if "NBPRJ_TEST_NBPATH" in os.environ:
        return _with_env(
            os.environ["NBPRJ_TEST_NBPATH"], env, return_env, default_env="test"
        )

    try:
        from IPython import get_ipython
    except ModuleNotFoundError:
        logger.warning("Can not import get_ipython.")
        return None

    main_module = sys.modules.get("__main__")
    if main_module is not None and hasattr(main_module, "__vsc_ipynb_file__"):
        return _with_env(
            main_module.__vsc_ipynb_file__,
            env,
            return_env,
            default_env="vs_code",
        )

    ipython_instance = get_ipython()
    if ipython_instance is None:
        logger.warning("The IPython instance is empty.")
        return None

    config = ipython_instance.config
    if "IPKernelApp" not in config:
        logger.warning("IPKernelApp is not in ipython_instance.config.")
        return None

    kernel_id = (
        config["IPKernelApp"]["connection_file"].partition("-")[2].split(".", -1)[0]
    )
    servers_nbapp, servers_juserv = _running_servers()
    server_exception = None
    for server in chain(servers_nbapp, servers_juserv):
        try:
            session = _query_server(server)
        except Exception as error:
            server_exception = error
            continue
        for notebook in session:
            if "kernel" not in notebook or "notebook" not in notebook:
                continue
            if notebook["kernel"].get("id", None) != kernel_id:
                continue
            for dir_key in _DIR_KEYS:
                if dir_key in server:
                    nb_path = PurePath(server[dir_key]) / notebook["notebook"]["path"]
                    default_env = "lab" if dir_key == "root_dir" else "notebook"
                    return _with_env(nb_path, env, return_env, default_env)

    nb_path = _lab_notebook_path()
    if nb_path is not None:
        return _with_env(nb_path, env, return_env, default_env="lab")

    # newer lab versions; stays unchanged after a file rename
    if "JPY_SESSION_NAME" in os.environ:
        return _with_env(
            PurePath(os.environ["JPY_SESSION_NAME"]),
            env,
            return_env,
            default_env="lab",
        )

    nb_path_psutil = _find_nb_path_via_parent_process()
    if nb_path_psutil is not None:
        logger.info("Detected path via psutil parent process inspection.")
        return _with_env(nb_path_psutil, env, return_env, default_env="nbconvert")

    if servers_nbapp == [] and servers_juserv == []:
        logger.warning("Can not find any servers running.")

    logger.warning(
        "Can not find the notebook in any server session or by using other methods."
    )
    if not _ipylab_is_installed():
        logger.warning(
            "Consider installing ipylab (pip install ipylab) if you use jupyter lab."
        )
    if server_exception is not None:
        raise server_exception
    return None
