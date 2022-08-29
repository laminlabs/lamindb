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

        if pipeline_run is not None:
            pipeline_run_id = pipeline_run.run_id
        else:
            pipeline_run_id = None

        with sqm.Session(engine) as session:
            result = session.get(db.schema.core.jupynb, (jupynb_id, jupynb_v))
            if result is None:
                session.add(
                    db.schema.core.jupynb(
                        id=jupynb_id,
                        v=jupynb_v,
                        name=jupynb_name,
                        user_id=settings.user.id,
                    )
                )
                dtransform_id = id.dtransform()
                session.add(
                    db.schema.core.dtransform(
                        id=dtransform_id,
                        jupynb_id=jupynb_id,
                        jupynb_v=jupynb_v,
                        pipeline_run_id=pipeline_run_id,
                    )
                )
                session.commit()
                logger.info(
                    f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                    f" user {settings.user.handle}."
                )
            else:
                dtransform = session.exec(
                    sqm.select(db.schema.core.dtransform).where(
                        db.schema.core.dtransform.jupynb_id == jupynb_id,
                        db.schema.core.dtransform.jupynb_v == jupynb_v,
                    )
                ).first()  # change to .one() as soon as dtransform ingestion bug fixed
                dtransform_id = dtransform.id

        with sqm.Session(engine) as session:
            if dobject_id is None:
                dobject_id = id.dobject()

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
            logger.success(f"Registered readout: {species.id}")

            return species.id

    @classmethod
    def features(
        cls,
        features_dict: dict,
        feature_entity: str,
        species: str,
        featureset_name: str = None,
        **kwargs,
    ):
        """Insert a featureset.

        Meanwhile inserting features and linking them to the featureset.
        """
        species_id = cls.species(common_name=species)

        # check if geneset exists
        if featureset_name is not None:
            query_featureset = getattr(db.do.query, "featureset")
            featureset_results = query_featureset(
                feature_entity=feature_entity,
                name=featureset_name,
            )
            if len(featureset_results) > 1:
                raise AssertionError(
                    f"Multiple entries are associated with {featureset_name}!"
                )
            elif len(featureset_results) == 1:
                logger.warning(f"Featureset {featureset_name} already exists!")
                return featureset_results[0].id

        # get the id field of feature entity
        feature_id = features_dict[next(iter(features_dict))].keys()[-1]
        query_feature = getattr(db.do.query, feature_entity)
        allfeatures = query_feature(species_id=species_id)
        # only ingest the new features but link all features to the featureset
        exist_feature_keys = set()
        exist_feature_ids = set()
        for feature in allfeatures:
            exist_feature_keys.add(feature.__getattribute__(feature_id))
            exist_feature_ids.add(feature.id)

        engine = settings.instance.db_engine()

        # add a featureset to the featureset table
        with sqm.Session(engine) as session:
            featureset = db.schema.bionty.featureset(
                feature_entity=feature_entity,
                name=featureset_name,
            )
            session.add(featureset)
            session.commit()
            session.refresh(featureset)

        # add features to the feature table
        with sqm.Session(engine) as session:
            features_ins = []
            for k, v in features_dict.items():
                if k in exist_feature_keys:
                    continue
                feature_schema = getattr(db.schema.bionty, feature_entity)
                feature = feature_schema(
                    **v,
                    species_id=species_id,
                )
                session.add(feature)
                features_ins.append(feature)
            session.commit()
            for feature in features_ins:
                session.refresh(feature)

        # insert ids into the link table
        feature_ids = [i.id for i in features_ins]
        feature_ids += exist_feature_ids
        with sqm.Session(engine) as session:
            for feature_id in feature_ids:
                featureset_link_module = getattr(
                    db.schema.bionty, f"featureset_{feature_entity}"
                )
                query_dict = {
                    "featureset_id": featureset.id,
                    f"{feature_entity}_id": feature_id,
                }
                link = featureset_link_module(**query_dict)
                session.add(link)
            session.commit()

        settings.instance._update_cloud_sqlite_file()

        return featureset.id

    @classmethod
    def readout(cls, efo_id: str):
        """Insert a row in the readout table."""
        assert sum(i.isdigit() for i in efo_id) == 7
        efo_id = efo_id.replace("_", ":")

        # check if entry already exists
        query_readout = getattr(db.do.query, "readout")
        readout_results = query_readout(efo_id=efo_id)
        if len(readout_results) > 1:
            raise AssertionError(f"Multiple entries are associated with {efo_id}!")
        elif len(readout_results) == 1:
            return readout_results[0].id
        else:
            engine = settings.instance.db_engine()

            from bioreadout import readout

            entry = readout(efo_id=efo_id)
            for k, v in entry.items():
                if isinstance(v, list):
                    entry[k] = ";".join(v)
            with sqm.Session(engine) as session:
                readout = db.schema.wetlab.readout(**entry)
                session.add(readout)
                session.commit()
                session.refresh(readout)
            logger.success(f"Registered readout: {readout.id}")

            return readout.id

    @classmethod
    def biometa(
        cls,
        dobject_id: str,
        biosample_id: int = None,
        readout_id: int = None,
        featureset_id: int = None,
    ):
        """Insert a row in the biometa table and link with a dobject."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            biometa = db.schema.wetlab.biometa(
                biosample_id=biosample_id,
                readout_id=readout_id,
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
    def pipeline(
        cls,
        id: str = None,
        v: str = None,
        name: str = None,
        reference: str = None,
    ):
        """Insert a new row in the pipeline table."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            pipeline = db.schema.core.pipeline(
                id=id, v=v, name=name, reference=reference
            )
            session.add(pipeline)
            session.commit()
            session.refresh(pipeline)

        settings.instance._update_cloud_sqlite_file()

        return pipeline.id

    @classmethod
    def pipeline_run(
        cls,
        id: str = None,
        name: str = None,
        pipeline_id: str = None,
        pipeline_v: str = None,
    ):
        """Insert a new row in the pipeline_run table."""
        engine = settings.instance.db_engine()

        with sqm.Session(engine) as session:
            pipeline_run = db.schema.core.pipeline_run(
                id=id, name=name, pipeline_id=pipeline_id, pipeline_v=pipeline_v
            )
            session.add(pipeline_run)
            session.commit()
            session.refresh(pipeline_run)

        settings.instance._update_cloud_sqlite_file()

        return pipeline_run.id
