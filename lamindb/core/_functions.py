import functools
import inspect
from contextvars import ContextVar
from datetime import datetime, timezone
from typing import Callable, Literal, ParamSpec, TypeVar

from lamin_utils import logger

from lamindb.base import deprecated

from ..models import Run, Transform
from ._context import (
    Context,
    detect_and_process_source_code_file,
    get_cli_args,
    serialize_params_to_json,
)
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
            # Get function metadata
            path, transform_type, reference, reference_type = (
                detect_and_process_source_code_file(
                    path=inspect.getsourcefile(func),
                    transform_type="function",
                )
            )

            initiated_by_run = get_current_tracked_run()
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

            # get the fully qualified module name, including submodules
            module_path = func.__module__.replace(".", "/")
            key = (
                module_path if module_path not in {"__main__", "__mp_main__"} else None
            )
            if path.exists():
                local_context = Context(uid=uid, path=path)
                local_context._create_or_load_transform(
                    description=None,
                    transform_type=transform_type,
                    transform_ref=reference,
                    transform_ref_type=reference_type,
                    key=key,
                )
                transform = local_context.transform
                transform._update_source_code_from_path(path)
                transform.save()
            else:
                if module_path in {"__main__", "__mp_main__"}:
                    qualified_name = f"{initiated_by_run.transform.key}"
                else:
                    qualified_name = f"{module_path}.py"
                transform = Transform(  # type: ignore
                    uid=uid,
                    key=qualified_name,
                    type="function",
                    source_code=inspect.getsource(func),
                ).save()

            run = Run(
                transform=transform,
                initiated_by_run=initiated_by_run,
                entrypoint=func.__qualname__,
            )  # type: ignore
            run.started_at = datetime.now(timezone.utc)
            run._status_code = -1  # started

            if is_flow:
                cli_args = get_cli_args()
                if cli_args:
                    logger.important(f"function invoked with: {cli_args}")
                    run.cli_args = cli_args

            # Bind arguments to get a mapping of parameter names to values
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            params = dict(bound_args.arguments)

            # Add parameters to the run
            run.params = serialize_params_to_json(params)
            run.save()

            # Set the run in context and execute function
            token = current_tracked_run.set(run)
            # If it's a flow, set the global run context as we do in `ln.track()`
            # Because we error above if a global run context already exists,
            # there is no danger of overwriting the global run context.
            if global_run in {"memorize", "clear"}:
                global_context._run = run
            try:
                result = func(*args, **kwargs)
                run.finished_at = datetime.now(timezone.utc)
                run._status_code = 0  # completed
                run.save()
                return result
            except Exception as e:
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
