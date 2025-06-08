"""helpers for providing lazy import of lamindb classes and funtions in lamindb.__init__."""

import sys
from types import new_class
from typing import Any, NoReturn

__all__ = [
    "new_metaclass_proxy",
    "new_instance_proxy",
    "new_callable_proxy",
]


def _raise_not_connected(objname: str) -> NoReturn:
    from lamindb_setup._check_setup import InstanceNotSetupError

    msg = f"To use lamindb.{objname}, you need to connect to an instance."
    raise InstanceNotSetupError(msg)


class _ProxyMeta(type):
    """MetaClass to enable lazy import of lamindb model classes."""

    def __call__(cls, *args, **kwargs):
        ln_cls = getattr(sys.modules["lamindb"], cls.__name__)
        if ln_cls is cls:
            _raise_not_connected(cls.__name__)
        return ln_cls(*args, **kwargs)

    def __subclasshook__(cls, __subclass):
        ln_cls = getattr(sys.modules["lamindb"], cls.__name__)
        if ln_cls is cls:
            _raise_not_connected(cls.__name__)
        return issubclass(__subclass, ln_cls)

    def __repr__(cls) -> str:
        ln_cls = getattr(sys.modules["lamindb"], cls.__name__)
        if ln_cls is cls:
            return f"<lamindb lazy_import proxy for {cls.__name__} (not connected)>"
        return repr(ln_cls)

    def __getattr__(cls, item):
        ln_cls = getattr(sys.modules["lamindb"], cls.__name__)
        if ln_cls is cls:
            _raise_not_connected(cls.__name__)
        return getattr(ln_cls, item)


class _Proxy:
    """Class proxy for lazily imported public instances in lamindb."""

    def __init__(self, name: str) -> None:
        self.__name = name

    def __repr__(self) -> str:
        obj = getattr(sys.modules["lamindb"], self.__name)
        if obj is None:
            return f"<lamindb lazy_import proxy for {self.__name} (not connected)>"
        return repr(obj)

    def __getattr__(self, item):
        obj = getattr(sys.modules["lamindb"], self.__name)
        if obj is None:
            _raise_not_connected(self.__name)
        return getattr(obj, item)

    def __setattr__(self, key, value):
        if key == "_Proxy__name":
            super().__setattr__(key, value)
        else:
            obj = getattr(sys.modules["lamindb"], self.__name)
            if obj is None:
                _raise_not_connected(self.__name)
            setattr(obj, key, value)


def new_metaclass_proxy(name: str) -> Any:
    """Create a metaclass proxy for top-level lamindb classes."""
    return new_class(name, (), {"metaclass": _ProxyMeta})


def new_instance_proxy(name: str) -> Any:
    """Create an instance proxy for top-level lamindb instances."""
    return _Proxy(name)


def new_callable_proxy(func_name: str) -> Any:
    """Create a callable proxy for top-level lamindb functions."""

    def func(*args, **kwargs):
        ln_func = getattr(sys.modules["lamindb"], func_name)
        if ln_func is func:
            _raise_not_connected(func_name)
        return ln_func(*args, **kwargs)

    func.__name__ = func_name
    return func
