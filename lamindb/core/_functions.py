import functools
import inspect
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Literal, ParamSpec, TypeVar

from lamindb.base import deprecated

from ..models import Run
from ._context import Context
from ._context import context as global_context

P = ParamSpec("P")
R = TypeVar("R")

# Create a context variable to store the current tracked run
current_tracked_run: ContextVar[Run | None] = ContextVar(
    "current_tracked_run", default=None
)


def get_current_tracked_run() -> Run | None:
    """Get the run object."""
    run = current_tracked_run.get()
    if run is None:
        run = global_context.run
    return run


def _create_tracked_decorator(
    uid: str | None = None,
    is_flow: bool = True,
    global_run: Literal["memorize", "clear", "none"] = "none",
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Internal helper to create tracked decorators.

    Args:
        uid: Persist the uid to identify this transform across renames.
        is_flow: Triggered through @ln.flow(), otherwise @ln.step().
    """

    def decorator_tracked(func: Callable[P, R]) -> Callable[P, R]:
        # Get the original signature
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper_tracked(*args: P.args, **kwargs: P.kwargs) -> R:
            if global_context.run is None:
                if not is_flow:
                    raise RuntimeError(
                        "Please track the global run context before using @ln.step(): ln.track() or @ln.flow()"
                    )
            else:
                if is_flow:
                    raise RuntimeError(
                        "Please clear the global run context before using @ln.flow(): no `ln.track()` or `@ln.flow(global_run='clear')`"
                    )
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            params = dict(bound_args.arguments)

            initiated_by_run = get_current_tracked_run()
            path_raw = inspect.getsourcefile(func)
            path = None
            # do not pass path when function is defined in an ipython cell
            if path_raw is not None and Path(path_raw).exists():
                path = Path(path_raw)
            module_path = func.__module__.replace(".", "/")
            key = (
                f"{module_path}.py"
                if module_path not in {"__main__", "__mp_main__"}
                else None
            )
            context = Context(uid=uid, path=path)
            context._track(
                path=path,
                key=key,
                source_code=inspect.getsource(func) if path is None else None,
                kind="function",
                entrypoint=func.__qualname__,
                params=params,
                new_run=True,
                initiated_by_run=initiated_by_run,
                stream_tracking=is_flow,
            )
            token = current_tracked_run.set(context.run)
            if global_run in {"memorize", "clear"}:
                global_context._run = context.run
            try:
                result = func(*args, **kwargs)
                context._finish()
                return result
            except Exception as e:
                run = context.run
                run.finished_at = datetime.now(timezone.utc)
                run._status_code = 1  # errored
                run.save()
                raise e
            finally:
                if (
                    global_run == "clear"
                    and global_context.run == current_tracked_run.get()
                ):
                    global_context._run = None
                current_tracked_run.reset(token)

        return wrapper_tracked

    return decorator_tracked


def flow(
    uid: str | None = None,
    global_run: Literal["memorize", "clear", "none"] = "clear",
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Use `@flow()` to track a function as a workflow.

    You will be able to see inputs, outputs, and parameters of the function in the data lineage graph.

    The decorator creates a :class:`~lamindb.Transform` object that maps onto the file in which the function is defined.
    The function maps onto an entrypoint of the `transform`.
    A function execution creates a :class:`~lamindb.Run` object that stores the function name in `run.entrypoint`.

    By default, like `ln.track()`, creates a global run context that can be accessed with `ln.context.run`.

    Args:
        uid: Persist the uid to identify a transform across renames.
        global_run: If `"clear"`, set the global run context `ln.context.run` and clear after the function completes.
            If `"memorize"`, set the global run context and do not clear after the function completes.
            Set this to `"none"` if you want to track concurrent executions of a `flow()` in the same Python process.

    Examples:

        To sync a workflow with a file in a git repo, see: :ref:`sync-code-with-git`.

        For an extensive guide, see: :ref:`manage-workflows`. Here follow some examples.

        .. literalinclude:: scripts/my_workflow.py
            :language: python
            :caption: my_workflow.py

        .. literalinclude:: scripts/my_workflow_with_step.py
            :language: python
            :caption: my_workflow_with_step.py

        .. literalinclude:: scripts/my_workflow_with_click.py
            :language: python
            :caption: my_workflow_with_click.py


    """
    return _create_tracked_decorator(uid=uid, is_flow=True, global_run=global_run)


def step(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Use `@step()` to track a function as a step.

    Behaves like :func:`~lamindb.flow()`, but acts as a step in a workflow and does not create a global run context.
    It errors if no initiating run (either global or local run context) exists.

    See :func:`~lamindb.flow()` for examples.

    Args:
        uid: Persist the uid to identify a transform across renames.
    """
    return _create_tracked_decorator(uid=uid, is_flow=False)


@deprecated("step")
def tracked(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    return step(uid)
