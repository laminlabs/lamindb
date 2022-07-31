import sqlmodel as sqm
from lamin_logger import logger
from lndb_schema_core import id
from lndb_setup import settings

import lamindb as db


class insert:
    """Insert data."""

    @classmethod
    def dobject_from_jupynb(
        cls,
        *,
        name: str,
        file_suffix: str = None,
        jupynb_id: str,
        jupynb_v: str,
        jupynb_name: str,
        dobject_id: str = None,
        dobject_v: str = "1",
    ):
        """Data object with its origin."""
        engine = settings.instance.db_engine()

        df_jupynb = db.do.load("jupynb")
        if jupynb_id not in df_jupynb.index:
            with sqm.Session(engine) as session:
                jupynb = db.schema.core.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    user_id=settings.user.user_id,
                )
                session.add(jupynb)
                session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.user_email} ({settings.user.user_id})."
            )

        with sqm.Session(engine) as session:
            dtransform_id = id.id_dtransform()
            dtransform = db.schema.core.dtransform(
                id=dtransform_id,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
            )
            session.add(dtransform)

            if dobject_id is None:
                dobject_id = id.id_dobject()

            dtransform_out = db.schema.core.dtransform_out(
                dtransform_id=dtransform_id,
                dobject_id=dobject_id,
                dobject_v=dobject_v,
            )
            session.add(dtransform_out)

            dobject = db.schema.core.dobject(
                id=dobject_id,
                v=dobject_v,
                name=name,
                dsource_id=dtransform_id,
                file_suffix=file_suffix,
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        settings.instance._update_cloud_sqlite_file()

        return dobject.id

    @classmethod
    def genes(
        cls,
        genes_dict: dict,
        geneset_name: str = None,
        species: str = None,
        **kwargs,
    ):
        """Insert a geneset.

        Mmeanwhile inserting genes and linking them to the geneset.
        """
        engine = settings.instance.db_engine()

        # add a geneset to the geneset table
        with sqm.Session(engine) as session:
            featureset = db.schema.bionty.featureset(
                feature_entity="gene",
                name=geneset_name,
            )
            session.add(featureset)
            session.commit()
            session.refresh(featureset)

        # add genes to the gene table
        with sqm.Session(engine) as session:
            genes_ins = []
            for _, v in genes_dict.items():
                gene = db.schema.bionty.gene(
                    **v,
                    species=species,
                )
                session.add(gene)
                genes_ins.append(gene)
            session.commit()
            for gene in genes_ins:
                session.refresh(gene)

        # insert ids into the link table
        with sqm.Session(engine) as session:
            for gene in genes_ins:
                link = db.schema.bionty.featureset_gene(
                    geneset_id=featureset.id,
                    gene_id=gene.id,
                )
                session.add(link)
            session.commit()

        settings.instance._update_cloud_sqlite_file()

        return featureset.id

    @classmethod
    def readout_type(cls, name: str, platform: str = None):
        """Insert a row in the readout table."""
        engine = settings.instance.db_engine()

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
        featureset_id: int = None,
    ):
        """Insert a row in the biometa table and link with a dobject."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            biometa = db.schema.wetlab.biometa(
                biosample_id=biosample_id,
                readout_type_id=readout_type_id,
                featureset_id=featureset_id,
            )
            session.add(biometa)
            session.commit()
            session.refresh(biometa)

        # also create an entry in the dobject_biometa table
        with sqm.Session(engine) as session:
            link = db.schema.wetlab.dobject_biometa(
                dobject_id=dobject_id, biometa_id=biometa.id
            )
            session.add(link)
            session.commit()
            session.refresh(link)

        settings.instance._update_cloud_sqlite_file()

        return biometa.id
