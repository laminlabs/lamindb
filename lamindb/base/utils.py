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
        # Copy metadata to the descriptor itself
        wraps(func)(self)

    def __get__(self, instance, owner):
        # Create a proper wrapper that preserves metadata
        if instance is None:

            @wraps(self.func)
            def wrapper(*args, **kwargs):
                return self.func(owner, *args, **kwargs)
        else:

            @wraps(self.func)
            def wrapper(*args, **kwargs):
                return self.func(instance, *args, **kwargs)

        return wrapper


class strict_classmethod:
    """Decorator for a classmethod that raises an error when called on an instance."""

    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __get__(self, instance, owner):
        if instance is not None:
            # Called on an instance - raise immediately or return cached error raiser
            raise TypeError(
                f"{owner.__name__}.{self.func.__name__}() is a class method and must be called on the {owner.__name__} class, not on a {owner.__name__} object"
            )

        # Called on the class - return bound method using MethodType
        if owner is None:
            owner = type(instance)
        return MethodType(self.func, owner)


__all__ = [
    "doc_args",
    "deprecated",
    "class_and_instance_method",
    "strict_classmethod",
]
