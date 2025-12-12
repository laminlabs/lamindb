"""Utilities.

.. autodecorator:: doc_args
.. autodecorator:: deprecated
.. autodecorator:: class_and_instance_method

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


__all__ = ["doc_args", "deprecated", "class_and_instance_method"]
