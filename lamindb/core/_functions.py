import functools
import inspect
from contextvars import ContextVar
from datetime import datetime, timezone
from typing import Callable, ParamSpec, TypeVar

from lamindb.base import deprecated

from ..models import Run, Transform
from ._context import context, serialize_params_to_json

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
        run = context.run
    return run


def _create_tracked_decorator(
    uid: str | None = None, raise_on_no_run: bool = True
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Internal helper to create tracked decorators.

    Args:
        uid: Persist the uid to identify this transform across renames.
        raise_on_no_run: If True, raise RuntimeError when no run context exists.
                         If False, skip tracking and execute function normally.
    """

    def decorator_tracked(func: Callable[P, R]) -> Callable[P, R]:
        # Get the original signature
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper_tracked(*args: P.args, **kwargs: P.kwargs) -> R:
            # Get function metadata
            source_code = inspect.getsource(func)

            initiated_by_run = get_current_tracked_run()
            if initiated_by_run is None:
                if context.run is None:
                    if raise_on_no_run:
                        raise RuntimeError(
                            "Please track the global run context before using @ln.step(): ln.track()"
                        )
                    else:
                        initiated_by_run = None
                else:
                    initiated_by_run = context.run
            # Get fully qualified function name
            module_name = func.__module__
            if module_name in {"__main__", "__mp_main__"}:
                qualified_name = (
                    f"{initiated_by_run.transform.key}/{func.__qualname__}.py"
                )
            else:
                qualified_name = f"{module_name}.{func.__qualname__}.py"

            # Create transform and run objects
            transform = Transform(  # type: ignore
                uid=uid,
                key=qualified_name,
                type="function",
                source_code=source_code,
            ).save()

            run = Run(transform=transform, initiated_by_run=initiated_by_run)  # type: ignore
            run.started_at = datetime.now(timezone.utc)
            run._status_code = -1  # started

            # Bind arguments to get a mapping of parameter names to values
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            params = dict(bound_args.arguments)

            # Add parameters to the run
            run.params = serialize_params_to_json(params)
            run.save()

            # Set the run in context and execute function
            token = current_tracked_run.set(run)
            try:
                result = func(*args, **kwargs)
                run.finished_at = datetime.now(timezone.utc)
                run._status_code = 0  # completed
                run.save()
                return result
            finally:
                current_tracked_run.reset(token)

        return wrapper_tracked

    return decorator_tracked


def step(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Track function runs.

    You will be able to see inputs, outputs, and parameters of the function in the data lineage graph.

    Guide: :doc:`/track`

    .. versionadded:: 1.1.0
        This is still in beta and will be refined in future releases.

    Args:
        uid: Persist the uid to identify this transform across renames.

    Example::

        import lamindb as ln

        @ln.step()
        def subset_dataframe(
            input_artifact_key: str,  # all arguments tracked as parameters of the function run
            output_artifact_key: str,
            subset_rows: int = 2,
            subset_cols: int = 2,
        ) -> None:
            artifact = ln.Artifact.get(key=input_artifact_key)
            df = artifact.load()  # auto-tracked as input
            new_df = df.iloc[:subset_rows, :subset_cols]
            ln.Artifact.from_dataframe(new_df, key=output_artifact_key).save()  # auto-tracked as output
    """
    return _create_tracked_decorator(uid=uid, raise_on_no_run=True)


def flow(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Track function runs without raising errors when no run context exists.

    Similar to :func:`step`, but gracefully skips tracking if no run context is available
    instead of raising a RuntimeError.

    Args:
        uid: Persist the uid to identify this transform across renames.

    Example::

        import lamindb as ln

        @ln.flow()
        def process_data(data: str) -> str:
            # This will track if run context exists, otherwise runs normally
            return data.upper()
    """
    return _create_tracked_decorator(uid=uid, raise_on_no_run=False)


@deprecated("step")
def tracked(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    return step(uid)
