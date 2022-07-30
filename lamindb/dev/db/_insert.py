from typing import Iterable, Union

import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings

import lamindb as db


class insert:
    """Insert data."""

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
        engine = settings.instance.db_engine()

        df_jupynb = db.do.load("jupynb")
        if jupynb_id not in df_jupynb.index:
            with sqm.Session(engine) as session:
                jupynb = db.schema.core.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    type=jupynb_type,
                    user_id=settings.user.user_id,
                )
                session.add(jupynb)
                session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.user_email} ({settings.user.user_id})."
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

        settings.instance._update_cloud_sqlite_file()

        return dobject.id

    @classmethod
    def readouts(
        cls,
        readouts: Iterable[str],
        readout_entity: str,
        readout_set_name: str = None,
        species: str = None,
        **kwargs,
    ):
        """Insert a geneset.

        Mmeanwhile inserting genes and linking them to the geneset.
        """
        engine = settings.instance.db_engine()
        if readout_entity == "gene":
            readoutset_schema = db.schema.bionty.geneset
            readout_schema = db.schema.bionty.gene
            link_schema = db.schema.bionty.geneset_gene
        elif readout_entity == "protein":
            readoutset_schema = db.schema.bionty.proteinset
            readout_schema = db.schema.bionty.protein
            link_schema = db.schema.bionty.proteinset_protein
        else:
            raise NotImplementedError

        # add a readout_set to the readout_set table
        with sqm.Session(engine) as session:
            readoutset = readoutset_schema(
                name=readout_set_name,
            )
            session.add(readoutset)
            session.commit()
            session.refresh(readoutset)

        # add readouts to the readout table
        with sqm.Session(engine) as session:
            readouts_ins = []
            for i in readouts:
                readout = readout_schema(
                    symbol=i,
                    species=species,
                    **kwargs,
                )
                session.add(readout)
                readouts_ins.append(readout)
            session.commit()
            for i in readouts_ins:
                session.refresh(i)

        # insert ids into the link table
        with sqm.Session(engine) as session:
            for readout in readouts_ins:
                link = link_schema(geneset_id=readoutset.id, gene_id=readout.id)
                session.add(link)
            session.commit()

        settings.instance._update_cloud_sqlite_file()

        return readoutset.id

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
        geneset_id: int = None,
        proteinset_id: int = None,
    ):
        """Insert a row in the biometa table and link with a dobject."""
        engine = settings.instance.db_engine()

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

        settings.instance._update_cloud_sqlite_file()

        return biometa.id

    @classmethod
    def biosample(cls, name: Union[str, Iterable[str]], species: str):
        """Insert entries in the biosample table."""
        engine = settings.instance.db_engine()

        names = [name] if isinstance(name, str) else name
        species_id = db.do.query.species(common_name=species).id
        with sqm.Session(engine) as session:
            biosamples: list = []
            for i in names:
                biosample = db.schema.biolab.biosample(name=i, species_id=species_id)
                session.add(biosample)
            session.commit()
            for i in biosamples:
                session.refresh(i)
        return [i.id for i in biosamples]

    @classmethod
    def species(
        cls, common_name: str, taxon_id: str, scientific_name: str, short_name: str
    ):
        """Insert a row in the species table."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            species = db.schema.bionty.species(
                common_name=common_name,
                taxon_id=taxon_id,
                scientific_name=scientific_name,
                short_name=short_name,
            )
            session.add(species)
            session.commit()
            session.refresh(species)

        return species.id
