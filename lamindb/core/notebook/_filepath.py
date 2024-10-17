import json
import os
from itertools import chain
from pathlib import PurePath
from urllib import request

from lamin_utils import logger

DIR_KEYS = ("notebook_dir", "root_dir")


def prepare_url(server: dict, query_str: str = ""):
    """Prepare url to query the jupyter server."""
    token = server["token"]
    if token:
        query_str = f"{query_str}?token={token}"
    url = f"{server['url']}api/sessions{query_str}"

    return url


def query_server(server: dict):
    """Query the jupyter server for sessions' info."""
    # based on https://github.com/msm1089/ipynbname
    try:
        url = prepare_url(server)
        with request.urlopen(url) as req:  # noqa: S310
            return json.loads(req.read())
    except Exception:
        CONN_ERROR = (
            "Unable to access server;\n"
            "querying requires either no security or token based security."
        )
        raise Exception(CONN_ERROR) from None


def running_servers():
    """Return the info about running jupyter servers."""
    nbapp_import = False

    try:
        from notebook.notebookapp import list_running_servers

        nbapp_import = True
        servers_nbapp = list_running_servers()
    except ModuleNotFoundError:
        servers_nbapp = []

    try:
        from jupyter_server.serverapp import list_running_servers

        servers_juserv = list_running_servers()
    except ModuleNotFoundError:
        servers_juserv = []

        if not nbapp_import:
            logger.warning(
                "It looks like you are running jupyter lab "
                "but don't have jupyter-server module installed. "
                "Please install it via pip install jupyter-server"
            )

    return servers_nbapp, servers_juserv


def notebook_path(return_env=False):
    """Return the path to the current notebook.

    Args:
        return_env: If `True`, return the environment of execution:
            `'lab'` for jupyter lab and `'notebook'` for jupyter notebook.
    """
    env = None
    if "NBPRJ_TEST_NBENV" in os.environ:
        env = os.environ["NBPRJ_TEST_NBENV"]

    if "NBPRJ_TEST_NBPATH" in os.environ:
        nb_path = os.environ["NBPRJ_TEST_NBPATH"]
        if return_env:
            return nb_path, "test" if env is None else env
        else:
            return nb_path

    try:
        from IPython import get_ipython
    except ModuleNotFoundError:
        logger.warning("Can not import get_ipython.")
        return None

    ipython_instance = get_ipython()

    # not in an ipython kernel
    if ipython_instance is None:
        logger.warning("The IPython instance is empty.")
        return None

    # vs code specific
    user_ns = ipython_instance.__dict__.get("user_ns", {})
    if "__vsc_ipynb_file__" in user_ns:
        nb_path = PurePath(user_ns["__vsc_ipynb_file__"])
        return (nb_path, "vs_code" if env is None else env) if return_env else nb_path

    config = ipython_instance.config
    # not in a jupyter notebook
    if "IPKernelApp" not in config:
        logger.warning("IPKernelApp is not in ipython_instance.config.")
        return None

    kernel_id = (
        config["IPKernelApp"]["connection_file"].partition("-")[2].split(".", -1)[0]
    )

    servers_nbapp, servers_juserv = running_servers()

    server_exception = None

    for server in chain(servers_nbapp, servers_juserv):
        try:
            session = query_server(server)
        except Exception as e:
            server_exception = e
            continue

        for notebook in session:
            if "kernel" not in notebook or "notebook" not in notebook:
                continue
            if notebook["kernel"].get("id", None) == kernel_id:
                for dir_key in DIR_KEYS:
                    if dir_key in server:
                        nb_path = (
                            PurePath(server[dir_key]) / notebook["notebook"]["path"]
                        )

                        if return_env:
                            if env is None:
                                rt_env = "lab" if dir_key == "root_dir" else "notebook"
                            else:
                                rt_env = env
                            return nb_path, rt_env
                        else:
                            return nb_path

    # for newer versions of lab, less safe as it stays the same after file rename
    if "JPY_SESSION_NAME" in os.environ:
        nb_path = PurePath(os.environ["JPY_SESSION_NAME"])
        return (nb_path, "lab" if env is None else env) if return_env else nb_path

    # no running servers
    if servers_nbapp == [] and servers_juserv == []:
        logger.warning("Can not find any servers running.")

    logger.warning(
        "Can not find the notebook in any server session or by using other methods."
    )

    if server_exception is not None:
        raise server_exception

    return None
