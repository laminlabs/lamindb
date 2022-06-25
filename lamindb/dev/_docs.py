from textwrap import dedent


def doc_args(*args):
    """Pass arguments to docstrings."""

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format(*args)
        return obj

    return dec
