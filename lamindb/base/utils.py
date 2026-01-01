"""Utilities.

.. autodecorator:: doc_args
.. autodecorator:: deprecated
.. autodecorator:: class_and_instance_method
.. autodecorator:: strict_classmethod

"""

from functools import wraps
from types import MethodType

from lamindb_setup.core import deprecated, doc_args


class class_and_instance_method:
    """Decorator to define a method that works both as class and instance method."""

    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __get__(self, instance, owner):
        if instance is None:
            # Called on the class
            return MethodType(self.func, owner)
        else:
            # Called on an instance
            return MethodType(self.func, instance)


class strict_classmethod:
    """Decorator for a classmethod that raises an error when called on an instance."""

    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __get__(self, instance, owner):
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
