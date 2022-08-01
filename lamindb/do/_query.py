from typing import Union

from lndb_setup import settings
from sqlmodel import Session, select
from sqlmodel.sql.expression import Select, SelectOfScalar

from .. import schema


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
        if (field == "cls") or (value is None):
            continue
        else:
            stmt = stmt.where(getattr(schema_module, field) == value)
    return stmt


class query:
    """Query literal (semantic) data."""

    @classmethod
    def id(cls, entity_name: str, id: Union[str, tuple]):
        """Query a single row by its id column with the primary key."""
        with Session(settings.instance.db_engine()) as session:
            for module in ["core", "wetlab", "bionty"]:
                schema_module = schema.__getattribute__(module)
                try:
                    return session.get(getattr(schema_module, entity_name), id)
                except AttributeError:
                    continue

    @classmethod
    def dobject(
        cls,
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
            featuresets = cls.featureset(gene=gene)
            biometas = []
            for i in featuresets:
                biometas += cls.biometa(featureset_id=i.id)
            dobjects = []
            for i in biometas:
                dobjects += cls.dobject_biometa(biometa_id=i.id)

        return [i for i in results if i.id in [j.dobject_id for j in dobjects]]

    @classmethod
    def user(cls, id: str = None, email: str = None, handle: str = None):
        """Query from the user table."""
        kwargs = locals()
        schema_module = schema.core.user
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    @classmethod
    def jupynb(
        cls, id: str = None, v: str = None, name: str = None, user_id: str = None
    ):
        """Query from the jupynb table."""
        kwargs = locals()
        schema_module = schema.core.user
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    @classmethod
    def readout_type(cls, id: int = None, name: str = None, platform: str = None):
        """Query from the readout_type table."""
        kwargs = locals()
        schema_module = schema.wetlab.readout_type
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    @classmethod
    def dobject_biometa(cls, dobject_id: str = None, biometa_id: str = None):
        """Query from the dobject_biometa table."""
        kwargs = locals()
        schema_module = schema.wetlab.dobject_biometa
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    @classmethod
    def biometa(
        cls,
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
            biometas = cls.dobject_biometa(dobject_id=dobject_id)
            if len(biometas) == 0:
                return biometas
            else:
                biometa_ids = [i.biometa_id for i in biometas]
                return [i for i in results if i.id in biometa_ids]
        else:
            return results

    @classmethod
    def featureset_gene(cls, featureset_id: int = None, gene_id: int = None):
        """Query the featureset_gene table."""
        kwargs = locals()
        schema_module = schema.bionty.featureset_gene
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        return _query_stmt(statement=stmt, results_type="all")

    @classmethod
    def featureset(
        cls,
        id: int = None,
        feature_entity: str = None,
        name: str = None,
        gene: str = None,
        protein: str = None,
    ):
        kwargs = locals()
        del kwargs["gene"]
        del kwargs["protein"]
        schema_module = schema.bionty.featureset
        stmt = _chain_select_stmt(kwargs=kwargs, schema_module=schema_module)
        results = _query_stmt(statement=stmt, results_type="all")
        if gene is not None:
            schema_module = schema.bionty.gene
            stmt = _chain_select_stmt(
                kwargs={"hgnc_symbol": gene}, schema_module=schema_module
            )
            gene_id = _query_stmt(statement=stmt, results_type="all")[0].id

            featuresets = cls.featureset_gene(gene_id=gene_id)
            featureset_ids = [i.featureset_id for i in featuresets]

        return [i for i in results if i.id in featureset_ids]
