import re

import pandas as pd
import sqlmodel as sqm
from lamin_logger import colors, logger
from lndb_setup import settings

from ...schema import core
from ...schema._table import Table
from ._select import select


def _camel_to_snake(string: str) -> str:
    """Convert CamelCase to snake_case."""

    def is_camel_case(s):
        return s != s.lower() and s != s.upper() and "_" not in s

    string = string.replace(" ", "_")
    if is_camel_case(string):
        return re.sub(r"(?<!^)(?=[A-Z])", "_", string).lower()
    return string.lower()


def dobject_from_dtransform(dobject: core.dobject, dtransform_id: str):
    storage = select.storage(  # type: ignore
        root=str(settings.instance.storage_root)
    ).one()

    dobject_id = insert.dobject(  # type: ignore
        id=dobject.id,
        name=dobject.name,
        dtransform_id=dtransform_id,
        suffix=dobject.suffix,
        storage_id=storage.id,
        size=dobject.size,
        checksum=dobject.checksum,
    )

    return dobject_id


def featureset_from_features(
    features_dict: dict,
    feature_entity: str,
    species: str,
    featureset_name: str = None,
):
    """Insert a featureset.

    Meanwhile inserting features and linking them to the featureset.
    """
    species = insert.species(common_name=species)  # type: ignore

    # check if geneset exists
    if featureset_name is not None:
        featureset_result = getattr(select, "featureset")(
            feature_entity=feature_entity,
            name=featureset_name,
        ).one_or_none()
        if featureset_result is not None:
            logger.warning(f"Featureset {featureset_name} already exists!")
            return featureset_result

    # get the id field of feature entity
    feature_id = features_dict[next(iter(features_dict))].keys()[-1]
    allfeatures = getattr(select, feature_entity)(species_id=species.id).all()  # type: ignore  # noqa
    # only ingest the new features but link all features to the featureset
    exist_feature_keys = set()
    exist_feature_ids = set()
    for feature in allfeatures:
        exist_feature_keys.add(feature.__getattribute__(feature_id))
        exist_feature_ids.add(feature.id)

    # add a featureset to the featureset table
    featureset = insert.featureset(  # type: ignore
        feature_entity=feature_entity, name=featureset_name
    )

    # add features to the feature table
    kwargs_list = []
    for k, v in features_dict.items():
        if k in exist_feature_keys:
            continue
        kwargs_list.append(v)
    added = InsertBase.insert_from_list(table_name=feature_entity, entries=kwargs_list)
    feature_ids = list(added.values()) + list(exist_feature_ids)
    for feature_id in feature_ids:
        kwargs = {
            "featureset_id": featureset.id,
            f"{feature_entity}_id": feature_id,
        }
        _ = getattr(insert, f"featureset_{feature_entity}")(**kwargs)

    return featureset


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
            results = cls.select(table_name=model.__name__, kwargs=kwargs)
            if len(results) >= 1:
                return "exists", results[0]
            elif len(results) > 1:
                return "exists", results

        with sqm.Session(settings.instance.db_engine()) as session:
            entry = model(**kwargs)
            session.add(entry)
            session.commit()
            session.refresh(entry)

        settings.instance._update_cloud_sqlite_file()

        return "inserted", entry

    @classmethod
    def select(cls, table_name, kwargs):
        return getattr(select, table_name)(**kwargs).all()

    @classmethod
    def is_unique(cls, model, column: str):
        return model.__table__.columns.get(column).unique

    @classmethod
    def insert_from_list(cls, table_name: str, entries):
        """Insert entries provided by a list of kwargs.

        Args:
            table_name: name of the table to insert
            entries: a list of table entries
        """
        if isinstance(entries[0], (dict, pd.Series)):
            model = Table.get_model(table_name)
            entries = [model(**d) for d in entries]

        added = {}
        with sqm.Session(settings.instance.db_engine()) as session:
            for i, entry in enumerate(entries):
                added[i] = entry
                session.add(added[i])
            session.commit()
            for i, _ in enumerate(entries):
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
        """Insert entries provided by a DataFrame.

        Raises a warning if not all columns in the schema table are passed.
        """
        mapper = {
            k: _camel_to_snake(k) for k in df.columns if k not in column_map.keys()
        }
        mapper.update(column_map)
        df = df.rename(columns=mapper).copy()

        # subset to columns that exist in the schema table
        fields = Table.get_fields(table_name)

        # if no mappable columns, raise an error
        fields_inters = set(fields).intersection(df.columns)
        if len(fields_inters) == 0:
            raise AssertionError(
                "No columns can be mapped between input DataFrame and table"
                f" {table_name}."
            )

        # if not all columns are populated, raise a warning
        fields_diff = set(fields).difference(df.columns)
        fields_diff.discard("id")
        fields_diff.discard("created_at")
        fields_diff.discard("updated_at")
        if len(fields_diff) > 0:
            logger.warning(f"The following columns are not populated: {fields_diff}")

        # insert entries into the table
        entries = df[list(fields_inters)].to_dict(orient="index")
        added = cls.insert_from_list(
            table_name=table_name, entries=list(entries.values())
        )
        entry_ids = {}
        for i, (idx, _) in enumerate(entries.items()):
            entry_ids[idx] = added[i]

        return entry_ids


def _create_insert_func(model):
    fields = Table.get_fields(model)
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
            status, entry = InsertBase.add(model=model, kwargs=kwargs, force=force)
        except AttributeError:
            status, entry = InsertBase.add(model=model, kwargs=kwargs, force=force)

        # no logging for these tables
        if status == "inserted" and name not in [
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
            entry_id = entry.id if hasattr(entry, "id") else entry
            logger.success(
                f"Inserted entry {colors.green(f'{entry_id}')} into"
                f" {colors.blue(f'{name}')}."
            )

        return entry

    insert_func.__name__ = name
    return insert_func


class insert:
    """Insert metadata.

    Guide: :doc:`/db/guide/insert-update-delete`.

    Example:

    >>> insert.pipeline(name="My pipeline A", v="1.3")

    Returns:
        id of the inserted entry
        None if entry already exists
    """

    pass


for model in Table.list_models():
    func = _create_insert_func(model=model)
    setattr(insert, model.__name__, staticmethod(func))

setattr(insert, "dobject_from_dtransform", staticmethod(dobject_from_dtransform))
setattr(insert, "featureset_from_features", staticmethod(featureset_from_features))
setattr(insert, "from_df", staticmethod(InsertBase.insert_from_df))
setattr(insert, "from_list", staticmethod(InsertBase.insert_from_list))
