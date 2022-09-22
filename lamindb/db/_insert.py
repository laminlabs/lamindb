import re
from typing import Iterable

import pandas as pd
import sqlmodel as sqm
from lamin_logger import colors, logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id
from sqlalchemy.orm.exc import NoResultFound

from .. import schema
from ..schema._table import Table
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


def features(
    features_dict: dict,
    feature_entity: str,
    species: str,
    featureset_name: str = None,
):
    """Insert a featureset.

    Meanwhile inserting features and linking them to the featureset.
    """
    species_id = getattr(insert, "species")(common_name=species)

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
    try:
        readout_id = getattr(query, "readout")(efo_id=efo_id).one()
    except NoResultFound:
        from bioreadout import readout

        entry = readout(efo_id=efo_id)
        readout_id = getattr(insert, "readout")(**entry)
        logger.success(
            f"Inserted entry {colors.green(f'{readout_id}')} into"
            f" {colors.blue('readout')}."
        )

    return readout_id


class FieldPopulator:
    @classmethod
    def species(cls, std_id_value: tuple) -> dict:
        from bionty import Species

        id_field, id_value = std_id_value
        ref_dict = Species(id=id_field).df.to_dict(orient="index")

        return ref_dict.get(id_value, {})

    @classmethod
    def readout(cls, std_id_value: tuple) -> dict:
        from bioreadout import readout

        id_field, id_value = std_id_value
        assert id_field == "efo_id"

        return readout(efo_id=id_field)


class InsertBase:
    @classmethod
    def add(cls, model, kwargs):
        with sqm.Session(settings.instance.db_engine()) as session:
            entry = model(**kwargs)
            session.add(entry)
            session.commit()
            session.refresh(entry)

        settings.instance._update_cloud_sqlite_file()

        return entry

    @classmethod
    def is_unique(cls, model, column):
        return model.__table__.columns.get(column).unique

    @classmethod
    def exists(cls, table_name, kwargs):
        results = getattr(query, table_name)(**kwargs).all()
        if len(results) == 0:
            return False
        return True

    @classmethod
    def insert_from_list(cls, entries: Iterable[dict], table_name: str):
        """Insert entries provided by a list of kwargs."""
        model = Table.get_model(table_name)
        added = {}
        with sqm.Session(settings.instance.db_engine()) as session:
            for i, kwargs in iter(entries):
                added[i] = model(**kwargs)
                session.add(added[i])
            session.commit()
            for i, kwargs in iter(entries):
                session.refresh(added[i])

        # fetch the ids
        if "id" in Table.get_pks(table_name):
            for k, v in added.items():
                added[k] = v.id
        else:
            for k, v in added.items():
                added[k] = v

        # returns {index : pk}
        return added

    @classmethod
    def insert_from_df(cls, df: pd.DataFrame, table_name: str, column_map: dict = {}):
        """Insert entries provided by a DataFrame."""
        mapper = {
            k: _camel_to_snake(k) for k in df.columns if k not in column_map.keys()
        }
        mapper.update(column_map)

        # subset to columns that exist in the schema table
        fields = Table.get_model(table_name).__fields__.keys()

        df = df.rename(columns=mapper).copy()
        df = df[df.columns.intersection(fields)]
        if df.shape[1] == 0:
            raise AssertionError(
                "No columns can be mapped between input DataFrame and table"
                f" {table_name}."
            )

        # insert entries into the table
        entries = df.to_dict(orient="index")
        entry_ids = {}
        for idx, entry in entries.items():
            entry_ids[idx] = getattr(insert, table_name)(**entry)

        return entry_ids


def _create_insert_func(table_name: str, model):
    def insert_func(cls, **kwargs):
        if InsertBase.exists(table_name=table_name, kwargs=kwargs):
            return
        try:
            reference = getattr(FieldPopulator, table_name)
            if len(kwargs) > 1:
                raise AssertionError(
                    "Please only provide a unique column id in the reference table."
                )
            std_id = next(iter(kwargs))
            std_value = kwargs[std_id]
            kwargs.update(reference(std_id_value=(std_id, std_value)))
            entry = InsertBase.add(model=model, kwargs=kwargs)
        except AttributeError:
            entry = InsertBase.add(model=model, kwargs=kwargs)

        pks = Table.get_pks(table_name)
        if "id" in pks:
            entry_id = entry.id
        else:
            entry_id = entry
        if table_name not in ["usage", "dobject"]:  # no logging
            logger.success(
                f"Inserted entry {colors.green(f'{entry_id}')} into"
                f" {colors.blue(f'{table_name}')}."
            )

        settings.instance._update_cloud_sqlite_file()

        return entry_id

    insert_func.__name__ = table_name
    return insert_func


class insert:
    """Insert an entry into a database table.

    Example:
    >>> insert.{entity}(id=1, name='new_experiment')
    """

    pass


for table_name, model in Table.all.items():
    func = _create_insert_func(table_name=table_name, model=model)
    setattr(insert, table_name, classmethod(func))

setattr(insert, "dobject_from_jupynb", dobject_from_jupynb)
setattr(insert, "features", features)
setattr(insert, "readout", readout)
setattr(insert, "from_df", InsertBase.insert_from_df)
setattr(insert, "from_list", InsertBase.insert_from_list)
