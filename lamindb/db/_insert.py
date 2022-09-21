import re

import pandas as pd
import sqlmodel as sqm
from lamin_logger import colors, logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

from .. import schema
from ..schema._schema import alltables
from ._query import query


def _camel_to_snake(string: str) -> str:
    """Convert CamelCase to snake_case."""

    def is_camel_case(s):
        return s != s.lower() and s != s.upper() and "_" not in s

    string = string.replace(" ", "_")
    if is_camel_case(string):
        return re.sub(r"(?<!^)(?=[A-Z])", "_", string).lower()
    return string.lower()


def dobject_from_jupynb(
    *,
    name: str,
    suffix: str = None,
    jupynb_id: str,
    jupynb_v: str,
    jupynb_name: str,
    dobject_id: str = None,
    dobject_v: str = "1",
    pipeline_run: BfxRun = None,
):
    """Data object from jupynb."""
    if pipeline_run is None:
        result = getattr(query, "jupynb")(id=jupynb_id, v=jupynb_v).all()
        if len(result) == 0:
            jupynb_id = getattr(insert, "jupynb")(
                id=jupynb_id,
                v=jupynb_v,
                name=jupynb_name,
                user_id=settings.user.id,
            )
            # dtransform entry
            dtransform_id = getattr(insert, "dtransform")(
                jupynb_id=jupynb_id, jupynb_v=jupynb_v
            )
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.handle}."
            )
        else:
            dtransform_id = (
                getattr(query, "dtransform")(jupynb_id=jupynb_id, jupynb_v=jupynb_v)
                .one()
                .id
            )
    else:
        result = getattr(query, "dtransform")(pipeline_run_id=pipeline_run.run_id).all()
        if len(result) == 0:
            dtransform_id = getattr(insert, "dtransform")(
                pipeline_run_id=pipeline_run.run_id
            )
        else:
            dtransform_id = result[0].id

    storage = getattr(query, "storage")(root=str(settings.instance.storage_dir)).first()
    if dobject_id is None:
        dobject_id = id.dobject()

    dobject_id = getattr(insert, "dobject")(
        id=dobject_id,
        v=dobject_v,
        name=name,
        dtransform_id=dtransform_id,
        suffix=suffix,
        storage_id=storage.id,
    )

    return dobject_id


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


def insert_from_df(df: pd.DataFrame, schema_table: str, column_map: dict = {}):
    """Insert entries provided by a DataFrame."""
    mapper = {k: _camel_to_snake(k) for k in df.columns if k not in column_map.keys()}
    mapper.update(column_map)

    # subset to columns that exist in the schema table
    table = alltables.get(schema_table)
    if table is None:
        raise AttributeError(f"Table {schema_table} is not found.")
    else:
        fields = table.__fields__.keys()

    df = df.rename(columns=mapper).copy()
    df = df[df.columns.intersection(fields)]
    if df.shape[1] == 0:
        raise AssertionError(
            "No columns can be mapped between input DataFrame and table"
            f" {schema_table}."
        )

    # insert entries into the table
    entries = df.to_dict(orient="index")
    entry_ids = {}
    for idx, entry in entries.items():
        entry_ids[idx] = getattr(insert, schema_table)(**entry)

    return entry_ids


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
        if name not in ["usage", "dobject"]:
            logger.success(
                f"Inserted entry {colors.green(f'{entry_id}')} into"
                f" {colors.blue(f'{name}')}."
            )

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
setattr(insert, "from_df", insert_from_df)
