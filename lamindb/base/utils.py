"""Utilities.

.. autodecorator:: doc_args
.. autodecorator:: deprecated
.. autodecorator:: class_and_instance_method
.. autofunction:: raise_error_if_called_on_object

"""

from functools import wraps

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


def raise_error_if_called_on_object(cls, method_name: str):
    """Raise an error if a classmethod is called on an object."""
    if hasattr(cls, "_state"):
        class_name = cls.__class__.__name__
        raise TypeError(
            f"{class_name}.{method_name}() is a class method and must be called on the {class_name} class, not on a {class_name} object"
        )


__all__ = [
    "doc_args",
    "deprecated",
    "class_and_instance_method",
    "raise_error_if_called_on_object",
]
