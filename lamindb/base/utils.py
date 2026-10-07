"""Utilities.

.. autodecorator:: doc_args
.. autodecorator:: deprecated
.. autodecorator:: class_and_instance_method
.. autodecorator:: strict_classmethod

"""

from collections.abc import Callable
from functools import wraps
from types import MethodType
from typing import (
    Any,
    Concatenate,
    Generic,
    NoReturn,
    ParamSpec,
    TypeVar,
    cast,
    overload,
)

from lamindb_setup.core import deprecated, doc_args

_T = TypeVar("_T")
_P = ParamSpec("_P")
_R = TypeVar("_R")


class class_and_instance_method(Generic[_T, _P, _R]):
    """Decorator to define a method that works both as class and instance method."""

    def __init__(self, func: Callable[Concatenate[_T, _P], _R]) -> None:
        self.func = func
        wraps(func)(cast(Any, self))

    @overload
    def __get__(self, instance: None, owner: type[_T]) -> Callable[_P, _R]: ...

    @overload
    def __get__(
        self, instance: _T, owner: type[_T] | None = None
    ) -> Callable[_P, _R]: ...

    def __get__(
        self, instance: _T | None, owner: type[_T] | None = None
    ) -> Callable[_P, _R]:
        if instance is None:
            # Called on the class
            if owner is None:
                raise TypeError("owner is required when accessing via class")
            return MethodType(self.func, owner)
        else:
            # Called on an instance
            return MethodType(self.func, instance)


class strict_classmethod(Generic[_T, _P, _R]):
    """Decorator for a classmethod that raises an error when called on an instance."""

    def __init__(self, func: Callable[Concatenate[type[_T], _P], _R]) -> None:
        self.func = func
        wraps(func)(cast(Any, self))

    @overload
    def __get__(self, instance: None, owner: type[_T]) -> Callable[_P, _R]: ...

    @overload
    def __get__(self, instance: _T, owner: type[_T] | None = None) -> NoReturn: ...

    def __get__(
        self, instance: _T | None, owner: type[_T] | None = None
    ) -> Callable[_P, _R]:
        if owner is None:
            raise TypeError("owner is required for descriptor access")
        if instance is not None:
            # Called on an instance - raise immediately
            raise TypeError(
                f"{owner.__name__}.{self.func.__name__}() is a class method and must be called on the {owner.__name__} class, not on a {owner.__name__} object"
            )

        # Called on the class - return bound method using MethodType
        return MethodType(self.func, owner)


def concrete_model(obj: Any) -> Any:
    """Return the model class whose table backs a record or model class.

    Resolves a Django proxy model to the concrete model it shares a table with, so
    a proxy of `Artifact` is dispatched like an `Artifact`. Returns any other class
    unchanged.

    The return type is `Any` because the result is the concrete model class of
    whatever was passed: a `Registry` stays a `Registry`, and an instance becomes
    that instance's model class (`objects`, `_meta`).
    """
    model = obj if isinstance(obj, type) else type(obj)
    concrete = getattr(getattr(model, "_meta", None), "concrete_model", None)
    return concrete or model


def get_registry_name(obj: Any) -> str:
    """Return the registry name of a record or model class, e.g. `"Artifact"`.

    Unlike `obj.__class__.__name__`, this is the same for a proxy model and the
    concrete model it proxies.
    """
    return concrete_model(obj).__name__


__all__ = [
    "doc_args",
    "deprecated",
    "class_and_instance_method",
    "strict_classmethod",
    "concrete_model",
    "get_registry_name",
]
