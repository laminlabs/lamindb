import sqlmodel as sqm
from lamin_logger import colors, logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

from .. import schema
from ..dev import track_usage
from ..schema._schema import alltables
from ._query import query


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
        result = session.get(schema.core.jupynb, (jupynb_id, jupynb_v))
        if result is None:
            session.add(
                schema.core.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    user_id=settings.user.id,
                )
            )
            dtransform_id = id.dtransform()
            session.add(
                schema.core.dtransform(
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
                sqm.select(schema.core.dtransform).where(
                    schema.core.dtransform.jupynb_id == jupynb_id,
                    schema.core.dtransform.jupynb_v == jupynb_v,
                )
            ).first()  # change to .one() as soon as dtransform ingestion bug fixed
            dtransform_id = dtransform.id

    with sqm.Session(engine) as session:
        if dobject_id is None:
            dobject_id = id.dobject()

        storage = session.exec(
            sqm.select(schema.core.storage).where(
                schema.core.storage.root == str(settings.instance.storage_dir)
            )
        ).first()
        assert storage

        dobject = schema.core.dobject(
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
    species_results = getattr(query, "species")(common_name=common_name).all()
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
            species = schema.bionty.species(**entry)
            session.add(species)
            session.commit()
            session.refresh(species)
        logger.success(
            f"Inserted entry {colors.green(f'{species.id}')} into"
            f" {colors.blue('species')}."
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
        featureset_results = getattr(query, "featureset")(
            feature_entity=feature_entity,
            name=featureset_name,
        ).all()
        if len(featureset_results) > 1:
            raise AssertionError(
                f"Multiple entries are associated with {featureset_name}!"
            )
        elif len(featureset_results) == 1:
            logger.warning(f"Featureset {featureset_name} already exists!")
            return featureset_results[0].id

    # get the id field of feature entity
    feature_id = features_dict[next(iter(features_dict))].keys()[-1]
    allfeatures = getattr(query, feature_entity)(species_id=species_id).all()
    # only ingest the new features but link all features to the featureset
    exist_feature_keys = set()
    exist_feature_ids = set()
    for feature in allfeatures:
        exist_feature_keys.add(feature.__getattribute__(feature_id))
        exist_feature_ids.add(feature.id)

    engine = settings.instance.db_engine()

    # add a featureset to the featureset table
    with sqm.Session(engine) as session:
        featureset = schema.bionty.featureset(
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
            feature_schema = getattr(schema.bionty, feature_entity)
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
                schema.bionty, f"featureset_{feature_entity}"
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
    readout_results = getattr(query, "readout")(efo_id=efo_id).all()
    if len(readout_results) > 1:
        raise AssertionError(f"Multiple entries are associated with {efo_id}!")
    elif len(readout_results) == 1:
        return readout_results[0].id
    else:
        from bioreadout import readout

        entry = readout(efo_id=efo_id)
        for k, v in entry.items():
            if isinstance(v, list):
                entry[k] = ";".join(v)
        with sqm.Session(settings.instance.db_engine()) as session:
            readout = schema.wetlab.readout(**entry)
            session.add(readout)
            session.commit()
            session.refresh(readout)
        logger.success(
            f"Inserted entry {colors.green(f'{readout.id}')} into"
            f" {colors.blue('readout')}."
        )

        settings.instance._update_cloud_sqlite_file()

        return readout.id


def _create_insert_func(name: str, schema_module):
    def insert_func(cls, **kwargs):
        with sqm.Session(settings.instance.db_engine()) as session:
            entry = schema_module(**kwargs)
            session.add(entry)
            session.commit()
            session.refresh(entry)
        try:
            entry_id = entry.id
        except AttributeError:
            entry_id = entry
        logger.success(
            f"Inserted entry {colors.green(f'{entry_id}')} into"
            f" {colors.blue(f'{name}')}."
        )
        if name == "dobject":
            track_usage(entry.id, entry.v, "insert")

        settings.instance._update_cloud_sqlite_file()

        return entry_id

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
