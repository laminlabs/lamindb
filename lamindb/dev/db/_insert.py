from typing import Iterable

import sqlmodel as sqm
from lamin_logger import logger
from lndb_cli import load_or_create_instance_settings, load_or_create_user_settings

import lamindb as db


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_email, user_id):
        df_user = db.do.load("user")
        if user_id not in df_user.index:
            user_id = insert.user(user_email, user_id)  # type: ignore
        return user_id


class insert:
    """Insert data."""

    @classmethod
    def schema_version(cls, version, user_id):
        """User."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()

        with sqm.Session(engine) as session:
            user = db.schema.core.schema_version(id=version, user_id=user_id)
            session.add(user)
            session.commit()

        settings._update_cloud_sqlite_file()

    @classmethod
    def user(cls, user_email, user_id):
        """User."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()

        with sqm.Session(engine) as session:
            user = db.schema.core.user(id=user_id, email=user_email)
            session.add(user)
            session.commit()
            session.refresh(user)

        settings._update_cloud_sqlite_file()

        return user.id

    @classmethod
    def dobject(
        cls,
        *,
        name: str,
        file_suffix: str = None,
        jupynb_id: str,
        jupynb_v: str,
        jupynb_name: str,
        jupynb_type: str,
        dobject_id: str = None,
        dobject_v: str = "1",
    ):
        """Data object with its origin."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()
        user_settings = load_or_create_user_settings()

        df_jupynb = db.do.load("jupynb")
        if jupynb_id not in df_jupynb.index:
            with sqm.Session(engine) as session:
                jupynb = db.schema.core.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    type=jupynb_type,
                    user_id=user_settings.user_id,
                )
                session.add(jupynb)
                session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {user_settings.user_email} ({user_settings.user_id})."
            )

        with sqm.Session(engine) as session:
            dobject = db.schema.core.dobject(
                id=dobject_id,
                v=dobject_v,
                name=name,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                file_suffix=file_suffix,
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        settings._update_cloud_sqlite_file()

        return dobject.id

    @classmethod
    def genes(
        cls,
        genes: Iterable[str],
        geneset_name: str = None,
        species: str = None,
        **kwargs,
    ):
        """Insert a geneset.

        Mmeanwhile inserting genes and linking them to the geneset.
        """
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()

        # add a geneset to the geneset table
        with sqm.Session(engine) as session:
            geneset = db.schema.bionty.geneset(
                name=geneset_name,
            )
            session.add(geneset)
            session.commit()
            session.refresh(geneset)

        # add genes to the gene table
        with sqm.Session(engine) as session:
            genes_ins = []
            for i in genes:
                gene = db.schema.bionty.gene(
                    symbol=i,
                    species=species,
                    **kwargs,
                )
                session.add(gene)
                genes_ins.append(gene)
            session.commit()
            for i in genes_ins:
                session.refresh(i)

        # insert ids into the link table
        with sqm.Session(engine) as session:
            for gene in genes_ins:
                link = db.schema.bionty.geneset_gene(
                    geneset_id=geneset.id,
                    gene_id=gene.id,
                )
                session.add(link)
            session.commit()

        settings._update_cloud_sqlite_file()

        return geneset.id

    @classmethod
    def readout_type(cls, name: str, platform: str = None):
        """Insert a row in the readout table."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()

        with sqm.Session(engine) as session:
            readout_type = db.schema.biolab.readout_type(name=name, platform=platform)
            session.add(readout_type)
            session.commit()
            session.refresh(readout_type)

        return readout_type.id

    @classmethod
    def biometa(
        cls,
        dobject_id: str,
        biosample_id: int = None,
        readout_type_id: int = None,
        geneset_id: int = None,
        proteinset_id: int = None,
    ):
        """Insert a row in the biometa table and link with a dobject."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()

        with sqm.Session(engine) as session:
            biometa = db.schema.biolab.biometa(
                biosample_id=biosample_id,
                readout_type_id=readout_type_id,
                geneset_id=geneset_id,
                proteinset_id=proteinset_id,
            )
            session.add(biometa)
            session.commit()
            session.refresh(biometa)

        # also create an entry in the dobject_biometa table
        with sqm.Session(engine) as session:
            link = db.schema.biolab.dobject_biometa(
                dobject_id=dobject_id, biometa_id=biometa.id
            )
            session.add(link)
            session.commit()
            session.refresh(link)

        return biometa.id
