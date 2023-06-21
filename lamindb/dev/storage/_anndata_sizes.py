import scipy.sparse as sparse
from pandas import DataFrame


def _size_val(val):
    if sparse.issparse(val):
        val_csr = sparse.csr_matrix(val)
        return val_csr.data.nbytes + val_csr.indptr.nbytes + val_csr.indices.nbytes
    else:
        return val.__sizeof__()


def _size_elem(elem):
    if hasattr(elem, "keys") and not isinstance(elem, DataFrame):
        return sum([_size_val(elem[k]) for k in elem.keys()])
    else:
        return _size_val(elem)


def _size_raw(raw):
    return _size_val(raw.X) + _size_val(raw.var) + _size_elem(raw.varm)


def size_adata(adata):
    total_size = adata.__sizeof__()

    raw = adata.raw
    if raw is not None:
        total_size += _size_raw(raw)

    return total_size
