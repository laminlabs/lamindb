import re
from typing import Iterable, Union

import pandas as pd
import sqlmodel as sqm
from lamin_logger import colors, logger
from lnbfx import BfxRun
from lndb_setup import settings
from lnschema_core import id

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


def dobject_from_dtransform(
    dobject_name: str,
    dtransform_id: str,
    dobject_suffix: Union[str, None],
    dobject_id: str,
    size: float,
):
    storage = getattr(query, "storage")(root=str(settings.instance.storage_dir)).first()
    if dobject_id is None:
        dobject_id = id.dobject()

    dobject_id = insert.dobject(  # type: ignore
        id=dobject_id,
        name=dobject_name,
        dtransform_id=dtransform_id,
        suffix=dobject_suffix,
        storage_id=storage.id,
        size=size,
    )

    return dobject_id


def dobject_from_pipeline(
    name: str,
    pipeline_run: BfxRun,
    suffix: Union[str, None],
    dobject_id: str,
    size: float,
):
    result = getattr(query, "dtransform")(
        pipeline_run_id=pipeline_run.run_id
    ).one_or_none()
    if result is None:
        dtransform_id = getattr(insert, "dtransform")(
            pipeline_run_id=pipeline_run.run_id
        )
    else:
        dtransform_id = result.id

    return dobject_from_dtransform(
        dobject_name=name,
        dobject_suffix=suffix,
        dobject_id=dobject_id,
        dtransform_id=dtransform_id,
        size=size,
    )


def dobject_from_jupynb(
    *,
    name: str,
    suffix: Union[str, None],
    jupynb_id: str,
    jupynb_v: str,
    jupynb_name: str,
    dobject_id: str,
    size: float,
):
    """Data object from jupynb."""
    result = getattr(query, "jupynb")(id=jupynb_id, v=jupynb_v).one_or_none()
    if result is None:
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

    return dobject_from_dtransform(
        dobject_name=name,
        dobject_suffix=suffix,
        dobject_id=dobject_id,
        dtransform_id=dtransform_id,
        size=size,
    )


def featureset_from_features(
    features_dict: dict,
    feature_entity: str,
    species: str,
    featureset_name: str = None,
):
    """Insert a featureset.

    Meanwhile inserting features and linking them to the featureset.
    """
    species_id = getattr(insert, "species")(common_name=species)
    if species_id is None:
        species_id = getattr(query, "species")(common_name=species).one().id

    # check if geneset exists
    if featureset_name is not None:
        featureset_result = getattr(query, "featureset")(
            feature_entity=feature_entity,
            name=featureset_name,
        ).one_or_none()
        if featureset_result is not None:
            logger.warning(f"Featureset {featureset_name} already exists!")
            return featureset_result.id

    # get the id field of feature entity
    feature_id = features_dict[next(iter(features_dict))].keys()[-1]
    allfeatures = getattr(query, feature_entity)(species_id=species_id).all()
    # only ingest the new features but link all features to the featureset
    exist_feature_keys = set()
    exist_feature_ids = set()
    for feature in allfeatures:
        exist_feature_keys.add(feature.__getattribute__(feature_id))
        exist_feature_ids.add(feature.id)

    # add a featureset to the featureset table
    featureset_id = getattr(insert, "featureset")(
        feature_entity=feature_entity, name=featureset_name
    )

    # add features to the feature table
    kwargs_list = []
    for k, v in features_dict.items():
        if k in exist_feature_keys:
            continue
        kwargs_list.append(v)
    added = InsertBase.insert_from_list(kwargs_list, feature_entity)
    feature_ids = list(added.values()) + list(exist_feature_ids)
    for feature_id in feature_ids:
        kwargs = {
            "featureset_id": featureset_id,
            f"{feature_entity}_id": feature_id,
        }
        _ = getattr(insert, f"featureset_{feature_entity}")(**kwargs)

    return featureset_id


class FieldPopulator:
    @classmethod
    def species(cls, std_id_value: tuple) -> dict:
        from bionty import Species

        id_field, id_value = std_id_value
        df = Species(id=id_field).df
        fields = Table.get_fields("species")
        df = df.loc[:, df.columns.intersection(fields)].copy()

        ref_dict = df.to_dict(orient="index")

        return ref_dict.get(id_value, {})

    @classmethod
    def readout(cls, std_id_value: tuple) -> dict:
        from bioreadout import readout

        id_field, id_value = std_id_value
        assert id_field == "efo_id"
        assert sum(i.isdigit() for i in id_value) == 7
        id_value = id_value.replace("_", ":")

        return readout(efo_id=id_value)


class InsertBase:
    @classmethod
    def add(cls, model, kwargs: dict, force=False):
        if not force:
            if cls.exists(table_name=model.__name__, kwargs=kwargs):
                return

        with sqm.Session(settings.instance.db_engine()) as session:
            entry = model(**kwargs)
            session.add(entry)
            session.commit()
            session.refresh(entry)

        settings.instance._update_cloud_sqlite_file()

        return entry

    @classmethod
    def exists(cls, table_name, kwargs):
        results = getattr(query, table_name)(**kwargs).all()
        if len(results) == 0:
            return False
        return True

    @classmethod
    def is_unique(cls, model, column: str):
        return model.__table__.columns.get(column).unique

    @classmethod
    def insert_from_list(cls, entries: Iterable[dict], table_name: str):
        """Insert entries provided by a list of kwargs."""
        model = Table.get_model(table_name)
        added = {}
        with sqm.Session(settings.instance.db_engine()) as session:
            for i, kwargs in enumerate(entries):
                added[i] = model(**kwargs)
                session.add(added[i])
            session.commit()
            for i, kwargs in enumerate(entries):
                session.refresh(added[i])

        # fetch the ids
        if "id" in Table.get_pks(table_name):
            for i, v in added.items():
                added[i] = v.id
        else:
            for i, v in added.items():
                added[i] = v

        settings.instance._update_cloud_sqlite_file()

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
        fields = Table.get_fields(table_name)

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


def _create_insert_func(model):
    fields = Table.get_fields(model)
    pks = Table.get_pks(model)
    name = model.__name__

    def insert_func(force=False, **kwargs):
        try:
            # insert with reference to populate the columns
            reference = getattr(FieldPopulator, name)
            # only allows one single kwarg
            if len(kwargs) > 1:
                raise AssertionError(
                    "Please only provide a unique column id in the reference table."
                )
            # check if the key is a unique column in the table
            std_id = next(iter(kwargs))
            if not InsertBase.is_unique(model, std_id):
                raise AssertionError("Please provide a unique column.")
            # populate other columns
            std_value = kwargs[std_id]
            kwargs_ = reference(std_id_value=(std_id, std_value))
            kwargs.update(**{k: v for k, v in kwargs_.items() if k in fields})
            entry = InsertBase.add(model=model, kwargs=kwargs, force=force)
        except AttributeError:
            entry = InsertBase.add(model=model, kwargs=kwargs, force=force)

        # returns None if an entry with the same kwargs already exists
        if entry is None:
            return

        # returns id or entry itself for link tables
        if "id" in pks:
            entry_id = entry.id
        else:
            entry_id = entry
        if name not in [
            "usage",
            "dobject",
            "gene",
            "protein",
            "cell_marker",
            "dobject_biometa",
            "dobject_bfxmeta",
            "featureset_gene",
            "featureset_protein",
            "featureset_cell_marker",
        ]:
            # no logging for these tables
            logger.success(
                f"Inserted entry {colors.green(f'{entry_id}')} into"
                f" {colors.blue(f'{name}')}."
            )

        return entry_id

    insert_func.__name__ = name
    return insert_func


class insert:
    """Insert an entry into a database table.

    Example:
    >>> insert.{entity}(id=1, name='new_experiment')

    Returns:
        id of the inserted entry
        None if entry already exists
    """

    pass


for model in Table.list_models():
    func = _create_insert_func(model=model)
    setattr(insert, model.__name__, staticmethod(func))

setattr(insert, "dobject_from_jupynb", staticmethod(dobject_from_jupynb))
setattr(insert, "dobject_from_pipeline", staticmethod(dobject_from_pipeline))
setattr(insert, "featureset_from_features", staticmethod(featureset_from_features))
setattr(insert, "from_df", staticmethod(InsertBase.insert_from_df))
setattr(insert, "from_list", staticmethod(InsertBase.insert_from_list))
