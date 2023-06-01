import operator
from typing import Union

# todo: add all operators
BINARY_OPS = [
    "__add__",
    "__mul__",
    "__or__",
    "__and__",
    "__eq__",
    "__ge__",
    "__gt__",
    "__le__",
    "__lt__",
    "__matmul__",
    "__pow__",
    "__truediv__",
    "__floordiv__",
]
UNARY_OPS = ["__abs__", "__neg__", "__invert__"]


class MetaCatchOperators(type):
    """Catch operators and return `LazyOperator`."""

    def __new__(meta, name, bases, dct):
        cls_obj = super().__new__(meta, name, bases, dct)

        for op in BINARY_OPS:
            op_func = getattr(operator, op)

            def left(self, right, op_func=op_func):
                return LazyOperator(self, right, op_func)

            setattr(cls_obj, op, left)

            def right(self, left, op_func=op_func):
                return LazyOperator(left, self, op_func)

            setattr(cls_obj, f"__r{op[2:]}", right)

        for op in UNARY_OPS:
            op_func = getattr(operator, op)

            def unary(self, op_func=op_func):
                return LazyOperator(self, None, op_func)

            setattr(cls_obj, op, unary)

        return cls_obj


class CatchAccess:
    def __getattr__(self, prop):
        """Catch a property and return `LazyProperty`."""
        return LazyProperty(self, prop)

    def __array_function__(self, func, types, args, kwargs):
        """Catch a numpy function and return `LazyNumpyFunc`."""
        return LazyNumpyFunc(func, args, kwargs)


class LazyOperator(CatchAccess, metaclass=MetaCatchOperators):
    def __init__(self, left, right, op):
        self._left = left
        self._right = right
        self._op = op

    def evaluate(self, **kwargs):
        _left_eval = self._left
        if hasattr(_left_eval, "evaluate"):
            _left_eval = _left_eval.evaluate(**kwargs)

        _right_eval = self._right
        if hasattr(_right_eval, "evaluate"):
            _right_eval = _right_eval.evaluate(**kwargs)

        if _right_eval is None:
            return self._op(_left_eval)
        else:
            return self._op(_left_eval, _right_eval)


class LazyProperty(CatchAccess, metaclass=MetaCatchOperators):
    def __init__(self, obj, name):
        self._name = name
        self._obj = obj
        self._call = False

        self._args = tuple()
        self._kwargs = dict()

    def evaluate(self, **kwargs):
        obj = self._obj
        if hasattr(self._obj, "evaluate"):
            obj = self._obj.evaluate(**kwargs)

        attr = getattr(obj, self._name)
        if self._call:
            result = attr(*self._args, **self._kwargs)
        else:
            result = attr
        return result

    def __call__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self._call = True

        return self


class LazyNumpyFunc(CatchAccess, metaclass=MetaCatchOperators):
    def __init__(self, func, args, kwargs):
        self._func = func
        self._args = args
        self._kwargs = kwargs

    def evaluate(self, **kwargs):
        args_eval = []
        for arg in self._args:
            if hasattr(arg, "evaluate"):
                arg = arg.evaluate(**kwargs)
            args_eval.append(arg)

        kwargs_eval = {}
        for key, val in self._kwargs.items():
            if hasattr(val, "evaluate"):
                val = val.evaluate(**kwargs)
            kwargs_eval[key] = val

        return self._func(*args_eval, **kwargs_eval)


class LazyField(CatchAccess, metaclass=MetaCatchOperators):
    def __init__(self, name, as_attr=True):
        self.name = name
        self._as_attr = as_attr

    def evaluate(self, obj):
        if self._as_attr:
            return getattr(obj, self.name)
        else:
            return obj[self.name]


class Lazy:
    def __getattr__(self, attr):
        """Get `LazyField` as an attribute."""
        return LazyField(attr, as_attr=True)

    def __getitem__(self, key):
        """Get `LazyField` as an item."""
        return LazyField(key, as_attr=False)


lazy = Lazy()
LazySelector = Union[LazyOperator, LazyProperty, LazyNumpyFunc]
