import functools
import inspect
from contextvars import ContextVar
from datetime import datetime, timezone
from typing import Callable, ParamSpec, TypeVar

from lamindb.base import deprecated

from ..models import Run, Transform
from ._context import (
    Context,
    detect_and_process_source_code_file,
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
    uid: str | None = None, is_flow: bool = True
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
                detect_and_process_source_code_file(path=inspect.getsourcefile(func))
            )

            initiated_by_run = get_current_tracked_run()
            if initiated_by_run is None:
                if global_context.run is None:
                    if not is_flow:
                        raise RuntimeError(
                            "Please track the global run context before using @ln.step(): ln.track()"
                        )
                    else:
                        initiated_by_run = None
                else:
                    initiated_by_run = global_context.run

            # Get fully qualified file name
            module_path = func.__module__.replace(
                ".", "/"
            )  # this is the fully qualified module name, including submodules
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
    return _create_tracked_decorator(uid=uid, is_flow=False)


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
    return _create_tracked_decorator(uid=uid, is_flow=True)


@deprecated("step")
def tracked(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    return step(uid)
