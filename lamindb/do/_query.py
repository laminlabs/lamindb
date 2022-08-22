from lndb_setup import settings
from sqlmodel import Session, select
from sqlmodel.sql.expression import Select, SelectOfScalar

from .. import schema
from ..schema._schema import alltables


def _query_stmt(statement, results_type="all"):
    with Session(settings.instance.db_engine()) as session:
        # Will remove after this is fixed:
        # https://github.com/tiangolo/sqlmodel/pull/234
        SelectOfScalar.inherit_cache = True  # type: ignore
        Select.inherit_cache = True  # type: ignore
        results = session.exec(statement).__getattribute__(results_type)()
    return results


def _chain_select_stmt(kwargs: dict, schema_module):
    stmt = select(schema_module)
    for field, value in kwargs.items():
        if (field in {"cls", "self"}) or (value is None):
            continue
        else:
            stmt = stmt.where(getattr(schema_module, field) == value)
    return stmt


def _create_query_func(name: str, schema_module):
    def query_func(cls, **kwargs):
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    query_func.__name__ = name
    return query_func


class query:
    """Query literal (semantic) data."""

    pass


def featureset(
    id: int = None,
    feature_entity: str = None,
    name: str = None,
    gene: str = None,
    protein: str = None,
):
    """Query the featureset table.

    Can also query a gene or a protein linked to featuresets.
    """
    kwargs = locals()
    del kwargs["gene"]
    del kwargs["protein"]
    schema_module = schema.bionty.featureset
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    if gene is not None:
        schema_module = schema.bionty.gene
        stmt = _chain_select_stmt(kwargs={"name": gene}, schema_module=schema_module)
        gene_id = _query_stmt(statement=stmt, results_type="all")[0].id

        featuresets = getattr(query, "featureset_gene")(gene_id=gene_id)
        featureset_ids = [i.featureset_id for i in featuresets]

        return [i for i in results if i.id in featureset_ids]
    else:
        return results


def biometa(
    id: int = None,
    biosample_id: int = None,
    readout_type_id: int = None,
    featureset_id: int = None,
    dobject_id: str = None,
):
    """Query from biometa.

    If dobject_id is provided, will search in the dobject_biometa first.
    """
    kwargs = locals()
    del kwargs["dobject_id"]
    schema_module = schema.wetlab.biometa
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    # dobject_id is given, will only return results associated with dobject_id
    if dobject_id is not None:
        biometas = getattr(query, "dobject_biometa")(dobject_id=dobject_id)
        if len(biometas) == 0:
            return biometas
        else:
            biometa_ids = [i.biometa_id for i in biometas]
            return [i for i in results if i.id in biometa_ids]
    else:
        return results


def dobject(
    id: str = None,
    v: str = None,
    name: str = None,
    file_suffix: str = None,
    dsource_id: str = None,
    gene: str = None,
):
    """Query from dobject."""
    kwargs = locals()
    del kwargs["gene"]
    schema_module = schema.core.dobject
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")

    if gene is not None:
        featuresets = getattr(query, "featureset")(gene=gene)
        biometas = []
        for i in featuresets:
            biometas += getattr(query, "biometa")(featureset_id=i.id)
        dobjects = []
        for i in biometas:
            dobjects += getattr(query, "dobject_biometa")(biometa_id=i.id)

        return [i for i in results if i.id in [j.dobject_id for j in dobjects]]
    else:
        return results


for name, schema_module in alltables.items():
    func = _create_query_func(name=name, schema_module=schema_module)
    setattr(query, name, classmethod(func))

setattr(query, "featureset", featureset)
setattr(query, "biometa", biometa)
setattr(query, "dobject", dobject)
