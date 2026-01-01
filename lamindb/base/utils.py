"""Utilities.

.. autodecorator:: doc_args
.. autodecorator:: deprecated
.. autodecorator:: class_and_instance_method
.. autodecorator:: raise_error_if_called_on_object

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


class raise_error_if_called_on_object:
    """Descriptor to raise an error if a classmethod is called on an instance."""

    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __get__(self, obj, objtype=None):
        print(f"DEBUG __get__ called for {self.func.__name__}")
        print(f"  obj: {obj}")
        print(f"  objtype: {objtype}")
        print(f"  obj is not None: {obj is not None}")

        if obj is not None:  # Called on an instance
            print("  -> Returning error_raiser for instance call")

            # Return a wrapper that will raise when called, not immediately
            def error_raiser(*args, **kwargs):
                class_name = objtype.__name__ if objtype else obj.__class__.__name__
                print(f"DEBUG error_raiser called: {class_name}.{self.func.__name__}")
                raise TypeError(
                    f"{class_name}.{self.func.__name__}() is a class method and must be called on the {class_name} class, not on a {class_name} object"
                )

            return error_raiser

        # Called on the class - return the bound classmethod
        print("  -> Returning bound classmethod for class call")
        return self.func.__get__(obj, objtype)

    def __call__(self, *args, **kwargs):
        print(f"DEBUG __call__ invoked for {self.func.__name__}")
        return self.func(*args, **kwargs)


__all__ = [
    "doc_args",
    "deprecated",
    "class_and_instance_method",
    "raise_error_if_called_on_object",
]
