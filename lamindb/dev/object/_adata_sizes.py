import scipy.sparse as sparse
from pandas import DataFrame


def _get_size_val(val):
    if sparse.issparse(val):
        val_csr = sparse.csr_matrix(val)
        return val_csr.data.nbytes + val_csr.indptr.nbytes + val_csr.indices.nbytes
    else:
        return val.__sizeof__()


def _get_size_elem(elem):
    if hasattr(elem, "keys") and not isinstance(elem, DataFrame):
        return sum([_get_size_val(elem[k]) for k in elem.keys()])
    else:
        return _get_size_val(elem)


def _get_size_raw(raw):
    return _get_size_val(raw.X) + _get_size_val(raw.var) + _get_size_elem(raw.varm)


def _get_size_adata(adata):
    total_size = adata.__sizeof__()

    raw = adata.raw
    if raw is not None:
        total_size += _get_size_raw(raw)

    return total_size
