import functools
import inspect
from contextvars import ContextVar
from datetime import datetime, timezone
from typing import Callable, ParamSpec, TypeVar

from .core._context import context
from .core._feature_manager import infer_feature_type_convert_json
from .models import Run, Transform

P = ParamSpec("P")
R = TypeVar("R")

# Create a context variable to store the current tracked run
current_tracked_run: ContextVar[Run | None] = ContextVar(
    "current_tracked_run", default=None
)


def get_current_tracked_run() -> Run | None:
    """Get the run object if we're inside a tracked function."""
    return current_tracked_run.get()


def tracked(uid: str | None = None) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Decorator that tracks function execution.

    Args:
        uid: Optional unique identifier for the transform
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
                    context.track()
                initiated_by_run = context.run

            # Get fully qualified function name
            module_name = func.__module__
            if module_name in {"__main__", "__mp_main__"}:
                qualified_name = f"{initiated_by_run.transform.key}/{func.__qualname__}"
            else:
                qualified_name = f"{module_name}.{func.__qualname__}"

            # Create transform and run objects
            transform = Transform(  # type: ignore
                uid=uid,
                key=qualified_name,
                type="function",
                source_code=source_code,
            ).save()

            run = Run(transform=transform, initiated_by_run=initiated_by_run)  # type: ignore
            run.started_at = datetime.now(timezone.utc)
            run.save()

            # Bind arguments to get a mapping of parameter names to values
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            params = dict(bound_args.arguments)

            # Remove the run parameter if it exists (we'll inject our own)
            params.pop("run", None)

            # Deal with non-trivial parameter values
            filtered_params = {}
            for key, value in params.items():
                dtype, _, _ = infer_feature_type_convert_json(
                    key, value, str_as_ulabel=False
                )
                if (dtype == "?" or dtype.startswith("cat")) and dtype != "cat ? str":
                    continue
                filtered_params[key] = value

            # Add parameters to the run
            run.params.add_values(filtered_params)

            # Set the run in context and execute function
            token = current_tracked_run.set(run)
            try:
                result = func(*args, **kwargs)
                run.finished_at = datetime.now(timezone.utc)
                run.save()
                return result
            finally:
                current_tracked_run.reset(token)

        return wrapper_tracked

    return decorator_tracked
