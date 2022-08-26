import sqlmodel as sqm
from lamin_logger import logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

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
        pipeline_run: BfxRun = None,
    ):
        """Data object with its origin."""
        engine = settings.instance.db_engine()

        df_jupynb = db.do.load.entity("jupynb")
        if jupynb_id not in df_jupynb.index:
            with sqm.Session(engine) as session:
                jupynb = db.schema.core.jupynb(
                    id=jupynb_id, v=jupynb_v, name=jupynb_name, user_id=settings.user.id
                )
                session.add(jupynb)
                session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.handle} ({settings.user.id})."
            )

        if pipeline_run is not None:
            pipeline_run_id = pipeline_run.run_id
        else:
            pipeline_run_id = None

        with sqm.Session(engine) as session:
            dtransform_id = id.id_dtransform()
            dtransform = db.schema.core.dtransform(
                id=dtransform_id,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                pipeline_run_id=pipeline_run_id,
            )
            session.add(dtransform)

            if dobject_id is None:
                dobject_id = id.id_dobject()

            storage = session.exec(
                sqm.select(db.schema.core.storage).where(
                    db.schema.core.storage.root == str(settings.instance.storage_dir)
                )
            ).first()
            assert storage

            dobject = db.schema.core.dobject(
                id=dobject_id,
                v=dobject_v,
                name=name,
                dtransform_id=dtransform_id,
                file_suffix=file_suffix,
                storage_id=storage.id,
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        settings.instance._update_cloud_sqlite_file()

        return dobject.id

    @classmethod
    def species(cls, common_name: str):
        """Insert a species."""
        query_species = getattr(db.do.query, "species")
        species_results = query_species(common_name=common_name)
        if len(species_results) > 1:
            raise AssertionError(f"Multiple entries are associated with {common_name}!")
        elif len(species_results) == 1:
            return species_results[0].id
        else:
            engine = settings.instance.db_engine()

            from bionty import Species

            entry = {"common_name": common_name}
            entry.update(Species().df.loc[common_name])
            with sqm.Session(engine) as session:
                species = db.schema.bionty.species(**entry)
                session.add(species)
                session.commit()
                session.refresh(species)
            logger.success(f"Registered readout_type: {species.id}")

            return species.id

    @classmethod
    def genes(
        cls,
        genes_dict: dict,
        species: str,
        geneset_name: str = None,
        **kwargs,
    ):
        """Insert a geneset.

        Meanwhile inserting genes and linking them to the geneset.
        """
        species_id = cls.species(common_name=species)

        # check if geneset exists
        if geneset_name is not None:
            query_featureset = getattr(db.do.query, "featureset")
            geneset_results = query_featureset(
                feature_entity="gene",
                name=geneset_name,
            )
            if len(geneset_results) > 1:
                raise AssertionError(
                    f"Multiple entries are associated with {geneset_name}!"
                )
            elif len(geneset_results) == 1:
                logger.warning(f"Geneset {geneset_name} already exists!")
                return geneset_results[0].id

        # get the id field
        gene_id = genes_dict[next(iter(genes_dict))].keys()[-1]
        query_gene = getattr(db.do.query, "gene")
        allgenes = query_gene(species=species_id)
        # only ingest the new genes but link all genes to the geneset
        exist_gene_keys = set()
        exist_gene_ids = set()
        for gene in allgenes:
            exist_gene_keys.add(gene.__getattribute__(gene_id))
            exist_gene_ids.add(gene.id)

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
            for k, v in genes_dict.items():
                if k in exist_gene_keys:
                    continue
                gene = db.schema.bionty.gene(
                    **v,
                    species=species_id,
                )
                session.add(gene)
                genes_ins.append(gene)
            session.commit()
            for gene in genes_ins:
                session.refresh(gene)

        # insert ids into the link table
        gene_ids = [i.id for i in genes_ins]
        gene_ids += exist_gene_ids
        with sqm.Session(engine) as session:
            for gene_id in gene_ids:
                link = db.schema.bionty.featureset_gene(
                    featureset_id=featureset.id,
                    gene_id=gene_id,
                )
                session.add(link)
            session.commit()

        settings.instance._update_cloud_sqlite_file()

        return featureset.id

    @classmethod
    def readout_type(cls, efo_id: str):
        """Insert a row in the readout table."""
        assert sum(i.isdigit() for i in efo_id) == 7
        efo_id = efo_id.replace("_", ":")

        # check if entry already exists
        query_readout = getattr(db.do.query, "readout_type")
        readout_results = query_readout(efo_id=efo_id)
        if len(readout_results) > 1:
            raise AssertionError(f"Multiple entries are associated with {efo_id}!")
        elif len(readout_results) == 1:
            return readout_results[0].id
        else:
            engine = settings.instance.db_engine()

            from bioreadout import readout_type

            entry = readout_type(efo_id=efo_id)
            for k, v in entry.items():
                if isinstance(v, list):
                    entry[k] = ";".join(v)
            with sqm.Session(engine) as session:
                readout_type = db.schema.wetlab.readout_type(**entry)
                session.add(readout_type)
                session.commit()
                session.refresh(readout_type)
            logger.success(f"Registered readout_type: {readout_type.id}")

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

    @classmethod
    def pipeline_run(cls, id: str = None):
        """Insert a new row in the pipeline_run table."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            pipeline_run = db.schema.core.pipeline_run(id=id)
            session.add(pipeline_run)
            session.commit()
            session.refresh(pipeline_run)

        settings.instance._update_cloud_sqlite_file()

        return pipeline_run.id
