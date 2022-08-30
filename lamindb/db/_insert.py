import sqlmodel as sqm
from lamin_logger import colors, logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

import lamindb as ln

from ..dev import track_usage
from ..schema._schema import alltables


def dobject_from_jupynb(
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
        result = session.get(ln.schema.core.jupynb, (jupynb_id, jupynb_v))
        if result is None:
            session.add(
                ln.schema.core.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    user_id=settings.user.id,
                )
            )
            dtransform_id = id.dtransform()
            session.add(
                ln.schema.core.dtransform(
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
                sqm.select(ln.schema.core.dtransform).where(
                    ln.schema.core.dtransform.jupynb_id == jupynb_id,
                    ln.schema.core.dtransform.jupynb_v == jupynb_v,
                )
            ).first()  # change to .one() as soon as dtransform ingestion bug fixed
            dtransform_id = dtransform.id

    with sqm.Session(engine) as session:
        if dobject_id is None:
            dobject_id = id.dobject()

        storage = session.exec(
            sqm.select(ln.schema.core.storage).where(
                ln.schema.core.storage.root == str(settings.instance.storage_dir)
            )
        ).first()
        assert storage

        dobject = ln.schema.core.dobject(
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


def insert_species(common_name: str):
    """Insert a species."""
    query_species = getattr(ln.db.query, "species")
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
            species = ln.schema.bionty.species(**entry)
            session.add(species)
            session.commit()
            session.refresh(species)
        logger.success(
            f"Inserted table {colors.blue('species')}: {colors.green(f'{species.id}')}"
        )

        return species.id


def features(
    features_dict: dict,
    feature_entity: str,
    species: str,
    featureset_name: str = None,
):
    """Insert a featureset.

    Meanwhile inserting features and linking them to the featureset.
    """
    species_id = insert_species(common_name=species)

    # check if geneset exists
    if featureset_name is not None:
        query_featureset = getattr(ln.db.query, "featureset")
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
    query_feature = getattr(ln.db.query, feature_entity)
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
        featureset = ln.schema.bionty.featureset(
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
            feature_schema = getattr(ln.schema.bionty, feature_entity)
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
                ln.schema.bionty, f"featureset_{feature_entity}"
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


def readout(efo_id: str):
    """Insert a row in the readout table."""
    assert sum(i.isdigit() for i in efo_id) == 7
    efo_id = efo_id.replace("_", ":")

    # check if entry already exists
    query_readout = getattr(ln.db.query, "readout")
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
            readout = ln.schema.wetlab.readout(**entry)
            session.add(readout)
            session.commit()
            session.refresh(readout)
        logger.success(
            f"Inserted table {colors.blue('readout')}: {colors.green(f'{readout.id}')}"
        )

        return readout.id


def biometa(
    dobject_id: str,
    biosample_id: int = None,
    readout_id: int = None,
    featureset_id: int = None,
):
    """Insert a row in the biometa table and link with a dobject."""
    engine = settings.instance.db_engine()

    with sqm.Session(engine) as session:
        biometa = ln.schema.wetlab.biometa(
            biosample_id=biosample_id,
            readout_id=readout_id,
            featureset_id=featureset_id,
        )
        session.add(biometa)
        session.commit()
        session.refresh(biometa)

    # also create an entry in the dobject_biometa table
    with sqm.Session(engine) as session:
        link = ln.schema.wetlab.dobject_biometa(
            dobject_id=dobject_id, biometa_id=biometa.id
        )
        session.add(link)
        session.commit()
        session.refresh(link)

    settings.instance._update_cloud_sqlite_file()

    return biometa.id


def _create_insert_func(name: str, schema_module):
    def insert_func(cls, **kwargs):
        with sqm.Session(settings.instance.db_engine()) as session:
            entry = schema_module(**kwargs)
            session.add(entry)
            session.commit()
            session.refresh(entry)
            logger.success(
                f"Inserted table {colors.blue(f'{name}')}:"
                f" {colors.green(f'{entry.id}')}"
            )
            if name == "dobject":
                track_usage(entry.id, entry.v, "insert")

        settings.instance._update_cloud_sqlite_file()

        return entry.id

    insert_func.__name__ = name
    return insert_func


class insert:
    """Insert an entry into a database table.

    Example:
    >>> insert.{entity}(id=1, name='new_experiment')
    """

    pass


for name, schema_module in alltables.items():
    func = _create_insert_func(name=name, schema_module=schema_module)
    setattr(insert, name, classmethod(func))

setattr(insert, "dobject_from_jupynb", dobject_from_jupynb)
setattr(insert, "species", insert_species)
setattr(insert, "features", features)
setattr(insert, "readout", readout)
setattr(insert, "biometa", biometa)
