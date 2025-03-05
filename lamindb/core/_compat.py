import importlib.util
from typing import Any, Callable, TypeVar

T = TypeVar("T")


def is_package_installed(package_name: str) -> bool:
    spec = importlib.util.find_spec(package_name)
    return spec is not None


def with_package(package_name: str, operation: Callable[[Any], T]) -> T:
    """Execute an operation that requires a specific package.

    Args:
        package_name: Package name (e.g., "mudata")
        operation: Function that takes the imported module and returns a result

    Examples:
        # For direct package functions
        result = with_package("mudata", lambda mod: mod.read_zarr(path))
    """
    try:
        module = importlib.import_module(package_name)
        return operation(module)
    except ImportError:
        raise ImportError(
            f"Package '{package_name}' is required but not installed. "
            f"Please install with: pip install {package_name}"
        ) from None


def with_package_obj(
    obj: Any, class_name: str, package_name: str, operation: Callable[[Any], T]
) -> tuple[bool, T | None]:
    """Handle operations on objects that require specific packages.

    Args:
        obj: The object to operate on
        class_name: Expected class name (e.g., "MuData")
        package_name: Package that provides the class (e.g., "mudata")
        operation: Function to call with the object if package is available.

    Examples:
        # For instance methods
        handled, res = apply_class_func(dmem, "MuData", "mudata",
                                      lambda obj: obj.write(filepath))
    """
    if obj.__class__.__name__ == class_name:
        try:
            importlib.import_module(package_name)
            result = operation(obj)
            return True, result
        except ImportError:
            raise ImportError(
                f"Object appears to be {class_name} but '{package_name}' package is not installed. "
                f"Please install with: pip install {package_name}"
            ) from None

    return False, None
