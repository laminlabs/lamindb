import bionty as bt
import bioreadout
import lnbfx


class lookup:
    """Lookup identifiers of knowledge tables."""

    pass


def _get_all_methods(lookup):
    return [i for i in dir(lookup) if not i.startswith("_")]


for module in [bt.lookup, bioreadout.lookup, lnbfx.lookup]:
    methods = _get_all_methods(module)
    for method in methods:
        model = getattr(module, method)
        setattr(lookup, method, staticmethod(model))
