import functools
import inspect
from datetime import datetime, timezone
from typing import Callable, ParamSpec, TypeVar

from ..models import Run, Transform
from ._context import context
from ._feature_manager import infer_feature_type_convert_json

P = ParamSpec("P")
R = TypeVar("R")


def tracked(
    uid: str | None = None, initiated_by_run: Run | None = None
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Decorator that tracks function execution in LaminDB and injects the run object.

    Args:
        uid: Optional unique identifier for the transform
        initiated_by_run: Optional parent run that initiated this function
    """

    def decorator_tracked(func: Callable[P, R]) -> Callable[P, R]:
        # Get the original signature
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper_tracked(*args: P.args, **kwargs: P.kwargs) -> R:
            nonlocal initiated_by_run
            # Get function metadata
            source_code = inspect.getsource(func)

            if initiated_by_run is None:
                assert context.run is not None  # noqa: S101
                initiated_by_run = context.run

            # Get fully qualified function name
            module_name = func.__module__
            if module_name == "__main__":
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
                dtype, _, _ = infer_feature_type_convert_json(key, value)
                if dtype == "?" or dtype.startswith("cat"):
                    continue
                filtered_params[key] = value

            # Add parameters to the run
            run.params.add_values(filtered_params)

            # Add the run to the kwargs
            kwargs["run"] = run

            # Call the original function with the injected run
            return func(*args, **kwargs)

        return wrapper_tracked

    return decorator_tracked
