import bionty as bt
import pandas as pd
from lndb_setup import settings
from sqlmodel import Session, select

from .. import schema
from ..dev import track_usage
from ..schema._schema import alltables


def _query_stmt(statement, results_type="all"):
    with Session(settings.instance.db_engine()) as session:
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


def _return_query_results_as_df(statement):
    df = pd.read_sql_query(statement, settings.instance.db_engine(future=False))
    if "id" in df.columns:
        if "v" in df.columns:
            df = df.set_index(["id", "v"])
        else:
            df = df.set_index("id")
    return df


def _create_query_func(name: str, schema_module):
    def query_func(cls, return_df=False, **kwargs):
        """Query metadata from tables."""
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        if return_df:
            results = _return_query_results_as_df(statement=stmt)
        else:
            results = _query_stmt(statement=stmt, results_type="all")

        # track usage for dobjects
        if name == "dobject":
            for result in results:
                track_usage(result.id, result.v, "query")
        return results

    query_func.__name__ = name
    return query_func


class query:
    """Query literal (semantic) data."""

    @classmethod
    def table_as_df(cls, entity_name) -> pd.DataFrame:
        """Load metadata table as DataFrame."""
        engine = settings.instance.db_engine()
        with engine.connect() as conn:
            df = pd.read_sql_table(entity_name, conn)
            if "id" in df.columns:
                if "v" in df.columns:
                    df = df.set_index(["id", "v"])
                else:
                    df = df.set_index("id")
        return df


# def featureset(id: int = None, feature_entity: str = None, name: str = None):
#     """Query the featureset table.

#     Can also query a gene or a protein linked to featuresets.
#     """
#     kwargs = locals()
#     del kwargs["gene"]
#     del kwargs["protein"]
#     del kwargs["cell_marker"]
#     schema_module = schema.bionty.featureset
#     stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
#     results = _query_stmt(statement=stmt, results_type="all")
#     if cell_marker is not None:
#         schema_module = schema.bionty.cell_marker
#         stmt = _chain_select_stmt(
#             kwargs={"name": cell_marker},  # TODO: remove hard code here
#             schema_module=schema_module,
#         )
#         feature_id = _query_stmt(statement=stmt, results_type="all")[0].id

#         featuresets = getattr(query, "featureset_cell_marker")(
#             cell_marker_id=feature_id
#         )
#         featureset_ids = [i.featureset_id for i in featuresets]

#         return [i for i in results if i.id in featureset_ids]
#     else:
#         return results


def _featureset_from_features(entity_name, **entity_kwargs):
    """Return featuresets by quering features."""
    schema_module = getattr(schema.bionty, entity_name)
    stmt = _chain_select_stmt(kwargs=entity_kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")
    if len(results) > 0:
        featureset_ids = []
        for feature in results:
            featuresets = getattr(query, f"featureset_{entity_name}")(
                **{f"{entity_name}_id": feature.id}
            )
            if len(featuresets) > 0:
                featureset_ids += [
                    featureset.featureset_id for featureset in featuresets
                ]
        return list(set(featureset_ids))
    else:
        return []


def biometa(
    id: int = None,
    biosample_id: int = None,
    readout_id: int = None,
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
    dtransform_id: str = None,
    file_suffix: str = None,
    storage_id: str = None,
    time_created=None,
    time_updated=None,
    entity_name: str = None,
    **entity_kwargs,
):
    """Query from dobject."""
    schema_module = schema.core.dobject
    kwargs = {
        k: v for k, v in locals().items() if k in schema.core.dobject.__fields__.keys()
    }
    stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
    results = _query_stmt(statement=stmt, results_type="all")

    if entity_name is not None:
        if getattr(bt.lookup.feature_model, entity_name) is not None:
            # query features
            featureset_ids = _featureset_from_features(
                entity_name=entity_name, **entity_kwargs
            )
            biometas = []
            for featureset_id in featureset_ids:
                biometas += getattr(query, "biometa")(featureset_id=featureset_id)
            dobjects = []
            for biometa in biometas:
                dobjects += getattr(query, "dobject_biometa")(biometa_id=biometa.id)
        else:
            # query obs metadata
            raise NotImplementedError

        results = [i for i in results if i.id in [j.dobject_id for j in dobjects]]

    if len(results) > 0:
        for result in results:
            track_usage(result.id, result.v, "query")

    return results


for name, schema_module in alltables.items():
    func = _create_query_func(name=name, schema_module=schema_module)
    setattr(query, name, classmethod(func))

setattr(query, "biometa", biometa)
setattr(query, "dobject", dobject)
