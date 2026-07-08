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


__all__ = [
    "doc_args",
    "deprecated",
    "class_and_instance_method",
    "strict_classmethod",
]
