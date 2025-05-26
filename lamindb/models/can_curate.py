from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Literal, Union

import numpy as np
import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from lamin_utils import colors, logger

from ..errors import ValidationError
from ._from_values import (
    _format_values,
    _from_values,
    get_organism_record_from_field,
)
from .sqlrecord import SQLRecord, get_name_field

if TYPE_CHECKING:
    from lamin_utils._inspect import InspectResult

    from lamindb.base.types import ListLike, StrField

    from .query_set import SQLRecordList


def _check_if_record_in_db(record: str | SQLRecord | None, using_key: str | None):
    """Check if the record is from the using_key DB."""
    if isinstance(record, SQLRecord):
        if using_key is not None and using_key != "default":
            if record._state.db != using_key:
                raise ValueError(
                    f"record must be a {record.__class__.__get_name_with_module__()} record from instance '{using_key}'!"
                )


def _concat_lists(values: ListLike | str) -> list[str]:
    """Concatenate a list of lists of strings into a single list."""
    if isinstance(values, str):
        values = [values]
    if isinstance(values, (list, pd.Series)) and len(values) > 0:
        first_item = values[0] if isinstance(values, list) else values.iloc[0]
        if isinstance(first_item, list):
            if isinstance(values, pd.Series):
                values = values.tolist()
            values = [
                v for sublist in values if isinstance(sublist, list) for v in sublist
            ]
    return values


def _inspect(
    cls,
    values: ListLike,
    field: StrField | None = None,
    *,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
    source: SQLRecord | None = None,
    from_source: bool = True,
    strict_source: bool = False,
) -> pd.DataFrame | dict[str, list[str]]:
    """{}"""  # noqa: D415
    from lamin_utils._inspect import inspect

    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)
    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    registry = queryset.model
    model_name = registry._meta.model.__name__
    if isinstance(source, SQLRecord):
        _check_if_record_in_db(source, queryset.db)
        # if strict_source mode, restrict the query to the passed ontology source
        # otherwise, inspect across records present in the DB from all ontology sources and no-source
        if strict_source:
            queryset = queryset.filter(source=source)
    organism_record = get_organism_record_from_field(
        getattr(registry, field_str), organism, values, queryset.db
    )
    _check_if_record_in_db(organism_record, queryset.db)

    # do not inspect synonyms if the field is not name field
    inspect_synonyms = True
    if hasattr(registry, "_name_field") and field_str != registry._name_field:
        inspect_synonyms = False

    # inspect in the DB
    result_db = inspect(
        df=_filter_queryset_with_organism(queryset=queryset, organism=organism_record),
        identifiers=values,
        field=field_str,
        mute=mute,
        inspect_synonyms=inspect_synonyms,
    )
    nonval = set(result_db.non_validated).difference(result_db.synonyms_mapper.keys())

    if from_source and len(nonval) > 0 and hasattr(registry, "source_id"):
        try:
            public_result = registry.public(
                organism=organism_record, source=source
            ).inspect(
                values=nonval,
                field=field_str,
                mute=True,
                inspect_synonyms=inspect_synonyms,
            )
            public_validated = public_result.validated
            public_mapper = public_result.synonyms_mapper
            hint = False
            if len(public_validated) > 0 and not mute:
                print_values = _format_values(public_validated)
                s = "" if len(public_validated) == 1 else "s"
                labels = colors.yellow(f"{len(public_validated)} {model_name} term{s}")
                logger.print(
                    f"   detected {labels} in public source for"
                    f" {colors.italic(field_str)}: {colors.yellow(print_values)}"
                )
                hint = True

            if len(public_mapper) > 0 and not mute:
                print_values = _format_values(list(public_mapper.keys()))
                s = "" if len(public_mapper) == 1 else "s"
                labels = colors.yellow(f"{len(public_mapper)} {model_name} term{s}")
                logger.print(
                    f"   detected {labels} in public source as {colors.italic(f'synonym{s}')}:"
                    f" {colors.yellow(print_values)}"
                )
                hint = True

            if hint:
                logger.print(
                    f"→  add records from public source to your {model_name} registry via"
                    f" {colors.italic('.from_values()')}"
                )

            nonval = [i for i in public_result.non_validated if i not in public_mapper]  # type: ignore
        # no public source is found
        except ValueError:
            logger.warning("no public source found, skipping source validation")

    if len(nonval) > 0 and not mute:
        print_values = _format_values(list(nonval))
        s = "" if len(nonval) == 1 else "s"
        labels = colors.red(f"{len(nonval)} term{s}")
        logger.print(f"   couldn't validate {labels}: {colors.red(print_values)}")
        logger.print(
            f"→  if you are sure, create new record{s} via"
            f" {colors.italic(f'{registry.__name__}()')} and save to your registry"
        )

    return result_db


def _validate(
    cls,
    values: ListLike,
    field: StrField | None = None,
    *,
    mute: bool = False,
    organism: str | SQLRecord | None = None,
    source: SQLRecord | None = None,
    strict_source: bool = False,
) -> np.ndarray:
    """{}"""  # noqa: D415
    from lamin_utils._inspect import validate

    return_str = True if isinstance(values, str) else False
    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)

    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    registry = queryset.model
    if isinstance(source, SQLRecord):
        _check_if_record_in_db(source, queryset.db)
        if strict_source:
            queryset = queryset.filter(source=source)

    organism_record = get_organism_record_from_field(
        getattr(registry, field_str), organism, values, queryset.db
    )
    _check_if_record_in_db(organism_record, queryset.db)
    field_values = pd.Series(
        _filter_queryset_with_organism(
            queryset=queryset,
            organism=organism_record,
            values_list_field=field_str,
        ),
        dtype="object",
    )
    if field_values.empty:
        if not mute:
            msg = (
                f"Your {cls.__name__} registry is empty, consider populating it first!"
            )
            if hasattr(cls, "source_id"):
                msg += "\n   → use `.import_source()` to import records from a source, e.g. a public ontology"
            logger.warning(msg)
        return np.array([False] * len(values))

    result = validate(
        identifiers=values,
        field_values=field_values,
        case_sensitive=True,
        mute=mute,
        field=field_str,
    )
    if return_str and len(result) == 1:
        return result[0]
    else:
        return result


def _standardize(
    cls,
    values: ListLike,
    field: StrField | None = None,
    *,
    return_field: str = None,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    mute: bool = False,
    source_aware: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    organism: str | SQLRecord | None = None,
    source: SQLRecord | None = None,
    strict_source: bool = False,
) -> list[str] | dict[str, str]:
    """{}"""  # noqa: D415
    from lamin_utils._standardize import standardize as map_synonyms

    return_str = True if isinstance(values, str) else False
    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)
    return_field_str = get_name_field(
        cls, field=field if return_field is None else return_field
    )
    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.objects.all()
    registry = queryset.model
    if isinstance(source, SQLRecord):
        _check_if_record_in_db(source, queryset.db)
        if strict_source:
            queryset = queryset.filter(source=source)
    organism_record = get_organism_record_from_field(
        getattr(registry, field_str), organism, values, queryset.db
    )
    _check_if_record_in_db(organism_record, queryset.db)

    # only perform synonym mapping if field is the name field
    if hasattr(registry, "_name_field") and field_str != registry._name_field:
        synonyms_field = None

    try:
        registry._meta.get_field(synonyms_field)
        fields = {
            field_name
            for field_name in [field_str, return_field_str, synonyms_field]
            if field_name is not None
        }
        df = _filter_queryset_with_organism(
            queryset=queryset,
            organism=organism_record,
            values_list_fields=list(fields),
        )
    except FieldDoesNotExist:
        df = pd.DataFrame()

    _kwargs = {
        "field": field_str,
        "return_field": return_field_str,
        "case_sensitive": case_sensitive,
        "keep": keep,
        "synonyms_field": synonyms_field,
    }
    # standardized names from the DB
    std_names_db = map_synonyms(
        df=df,
        identifiers=values,
        return_mapper=return_mapper,
        mute=mute,
        **_kwargs,
    )

    def _return(result: list, mapper: dict):
        if return_mapper:
            return mapper
        else:
            if return_str and len(result) == 1:
                return result[0]
            return result

    # map synonyms in public source
    if hasattr(registry, "source_id") and source_aware:
        mapper = {}
        if return_mapper:
            mapper = std_names_db
            std_names_db = map_synonyms(
                df=df, identifiers=values, return_mapper=False, mute=True, **_kwargs
            )

        val_res = registry.validate(
            std_names_db, field=field, mute=True, organism=organism_record
        )
        if all(val_res):
            return _return(result=std_names_db, mapper=mapper)

        nonval = np.array(std_names_db)[~val_res]
        std_names_bt_mapper = registry.public(organism=organism_record).standardize(
            nonval, return_mapper=True, mute=True, **_kwargs
        )

        if len(std_names_bt_mapper) > 0 and not mute:
            s = "" if len(std_names_bt_mapper) == 1 else "s"
            field_print = "synonym" if field_str == return_field_str else field_str

            reduced_mapped_keys_str = f"{list(std_names_bt_mapper.keys())[:10] + ['...'] if len(std_names_bt_mapper) > 10 else list(std_names_bt_mapper.keys())}"
            truncated_note = (
                " (output truncated)" if len(std_names_bt_mapper) > 10 else ""
            )

            warn_msg = (
                f"found {len(std_names_bt_mapper)} {field_print}{s} in public source{truncated_note}:"
                f" {reduced_mapped_keys_str}\n"
                f"  please add corresponding {registry._meta.model.__name__} records via{truncated_note}:"
                f" `.from_values({reduced_mapped_keys_str})`"
            )

            logger.warning(warn_msg)

        mapper.update(std_names_bt_mapper)
        if isinstance(std_names_db, pd.CategoricalDtype):
            result = std_names_db.cat.rename_categories(std_names_bt_mapper).tolist()
        else:
            result = pd.Series(std_names_db).replace(std_names_bt_mapper).tolist()
        return _return(result=result, mapper=mapper)

    else:
        return _return(result=std_names_db, mapper=std_names_db)


def _add_or_remove_synonyms(
    synonym: str | ListLike,
    record: CanCurate,
    action: Literal["add", "remove"],
    force: bool = False,
    save: bool | None = None,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: set[str], record: CanCurate):
        """Errors if input synonym is associated with other records in the DB."""
        import pandas as pd
        from IPython.display import display

        syns_all = (
            record.__class__.objects.exclude(synonyms="").exclude(synonyms=None).all()  # type: ignore
        )
        if len(syns_all) == 0:
            return
        df = pd.DataFrame(syns_all.values())
        df["synonyms"] = df["synonyms"].str.split("|")
        df = df.explode("synonyms")
        matches_df = df[(df["synonyms"].isin(synonyms)) & (df["id"] != record.id)]  # type: ignore
        if matches_df.shape[0] > 0:
            records_df = pd.DataFrame(syns_all.filter(id__in=matches_df["id"]).values())
            logger.error(
                f"input synonyms {matches_df['synonyms'].unique()} already associated"
                " with the following records:\n"
            )
            display(records_df)
            raise ValidationError(
                f"you are trying to assign a synonym to record: {record}\n"
                "    → consider removing the synonym from existing records or using a different synonym."
            )

    # passed synonyms
    # nothing happens when passing an empty string or list
    if isinstance(synonym, str):
        if len(synonym) == 0:
            return
        syn_new_set = {synonym}
    else:
        if synonym == [""]:
            return
        syn_new_set = set(synonym)
    # nothing happens when passing an empty string or list
    if len(syn_new_set) == 0:
        return
    # because we use | as the separator
    if any("|" in i for i in syn_new_set):
        raise ValidationError("a synonym can't contain '|'!")

    # existing synonyms
    syns_exist = record.synonyms  # type: ignore
    if syns_exist is None or len(syns_exist) == 0:
        syns_exist_set = set()
    else:
        syns_exist_set = set(syns_exist.split("|"))

    if action == "add":
        if not force:
            check_synonyms_in_all_records(syn_new_set, record)
        syns_exist_set.update(syn_new_set)
    elif action == "remove":
        syns_exist_set = syns_exist_set.difference(syn_new_set)

    if len(syns_exist_set) == 0:
        syns_str = None
    else:
        syns_str = "|".join(syns_exist_set)

    record.synonyms = syns_str  # type: ignore

    if save is None:
        # if record is already in DB, save the changes to DB
        save = not record._state.adding  # type: ignore
    if save:
        record.save()  # type: ignore


def _check_synonyms_field_exist(record: CanCurate):
    """Check if synonyms field exists."""
    if not hasattr(record, "synonyms"):
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        ) from None


def _filter_queryset_with_organism(
    queryset: QuerySet,
    organism: SQLRecord | None = None,
    values_list_field: str | None = None,
    values_list_fields: list[str] | None = None,
):
    """Filter a queryset based on organism."""
    import pandas as pd

    if organism is not None:
        queryset = queryset.filter(organism=organism)

    # values_list_field/s for better performance
    if values_list_field is None:
        if values_list_fields:
            return pd.DataFrame.from_records(
                queryset.values_list(*values_list_fields), columns=values_list_fields
            )
        return pd.DataFrame.from_records(queryset.values())
    else:
        return queryset.values_list(values_list_field, flat=True)


class CanCurate:
    """Base class providing :class:`~lamindb.models.SQLRecord`-based validation."""

    @classmethod
    def inspect(
        cls,
        values: ListLike,
        field: StrField | None = None,
        *,
        mute: bool = False,
        organism: Union[str, SQLRecord, None] = None,
        source: SQLRecord | None = None,
        from_source: bool = True,
        strict_source: bool = False,
    ) -> InspectResult:
        """Inspect if values are mappable to a field.

        Being mappable means that an exact match exists.

        Args:
            values: Values that will be checked against the field.
            field: The field of values. Examples are `'ontology_id'` to map
                against the source ID or `'name'` to map against the ontologies
                field names.
            mute: Whether to mute logging.
            organism: An Organism name or record.
            source: A `bionty.Source` record that specifies the version to inspect against.
            strict_source: Determines the validation behavior against records in the registry.
                - If `False`, validation will include all records in the registry, ignoring the specified source.
                - If `True`, validation will only include records in the registry  that are linked to the specified source.
                Note: this parameter won't affect validation against public sources.

        See Also:
            :meth:`~lamindb.models.CanCurate.validate`

        Example::

            import bionty as bt

            # save some gene records
            bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol", organism="human").save()

            # inspect gene symbols
            gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            result = bt.Gene.inspect(gene_symbols, field=bt.Gene.symbol, organism="human")
            assert result.validated == ["A1CF", "A1BG"]
            assert result.non_validated == ["FANCD1", "FANCD20"]
        """
        return _inspect(
            cls=cls,
            values=values,
            field=field,
            mute=mute,
            strict_source=strict_source,
            organism=organism,
            source=source,
            from_source=from_source,
        )

    @classmethod
    def validate(
        cls,
        values: ListLike,
        field: StrField | None = None,
        *,
        mute: bool = False,
        organism: Union[str, SQLRecord, None] = None,
        source: SQLRecord | None = None,
        strict_source: bool = False,
    ) -> np.ndarray:
        """Validate values against existing values of a string field.

        Note this is strict_source validation, only asserts exact matches.

        Args:
            values: Values that will be validated against the field.
            field: The field of values.
                    Examples are `'ontology_id'` to map against the source ID
                    or `'name'` to map against the ontologies field names.
            mute: Whether to mute logging.
            organism: An Organism name or record.
            source: A `bionty.Source` record that specifies the version to validate against.
            strict_source: Determines the validation behavior against records in the registry.
                - If `False`, validation will include all records in the registry, ignoring the specified source.
                - If `True`, validation will only include records in the registry  that are linked to the specified source.
                Note: this parameter won't affect validation against public sources.

        Returns:
            A vector of booleans indicating if an element is validated.

        See Also:
            :meth:`~lamindb.models.CanCurate.inspect`

        Example::

            import bionty as bt

            bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol", organism="human").save()

            gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            bt.Gene.validate(gene_symbols, field=bt.Gene.symbol, organism="human")
            #> array([ True,  True, False, False])
        """
        return _validate(
            cls=cls,
            values=values,
            field=field,
            mute=mute,
            strict_source=strict_source,
            organism=organism,
            source=source,
        )

    @classmethod
    def from_values(
        cls,
        values: ListLike,
        field: StrField | None = None,
        create: bool = False,
        organism: Union[SQLRecord, str, None] = None,
        source: SQLRecord | None = None,
        mute: bool = False,
    ) -> SQLRecordList:
        """Bulk create validated records by parsing values for an identifier such as a name or an id).

        Args:
            values: A list of values for an identifier, e.g.
                `["name1", "name2"]`.
            field: A `SQLRecord` field to look up, e.g., `bt.CellMarker.name`.
            create: Whether to create records if they don't exist.
            organism: A `bionty.Organism` name or record.
            source: A `bionty.Source` record to validate against to create records for.
            mute: Whether to mute logging.

        Returns:
            A list of validated records. For bionty registries. Also returns knowledge-coupled records.

        Notes:
            For more info, see tutorial: :doc:`docs:bio-registries`.

        Example::

            import bionty as bt

            # Bulk create from non-validated values will log warnings & returns empty list
            ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"])
            assert len(ulabels) == 0

            # Bulk create records from validated values returns the corresponding existing records
            ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"], create=True).save()
            assert len(ulabels) == 3

            # Bulk create records from public reference
            bt.CellType.from_values(["T cell", "B cell"]).save()
        """
        return _from_values(
            iterable=values,
            field=getattr(cls, get_name_field(cls, field=field)),
            create=create,
            organism=organism,
            source=source,
            mute=mute,
        )

    @classmethod
    def standardize(
        cls,
        values: Iterable,
        field: StrField | None = None,
        *,
        return_field: StrField | None = None,
        return_mapper: bool = False,
        case_sensitive: bool = False,
        mute: bool = False,
        source_aware: bool = True,
        keep: Literal["first", "last", False] = "first",
        synonyms_field: str = "synonyms",
        organism: Union[str, SQLRecord, None] = None,
        source: SQLRecord | None = None,
        strict_source: bool = False,
    ) -> list[str] | dict[str, str]:
        """Maps input synonyms to standardized names.

        Args:
            values: Identifiers that will be standardized.
            field: The field representing the standardized names.
            return_field: The field to return. Defaults to field.
            return_mapper: If `True`, returns `{input_value: standardized_name}`.
            case_sensitive: Whether the mapping is case sensitive.
            mute: Whether to mute logging.
            source_aware: Whether to standardize from public source. Defaults to `True` for BioRecord registries.
            keep: When a synonym maps to multiple names, determines which duplicates to mark as `pd.DataFrame.duplicated`:
                - `"first"`: returns the first mapped standardized name
                - `"last"`: returns the last mapped standardized name
                - `False`: returns all mapped standardized name.

                When `keep` is `False`, the returned list of standardized names will contain nested lists in case of duplicates.

                When a field is converted into return_field, keep marks which matches to keep when multiple return_field values map to the same field value.
            synonyms_field: A field containing the concatenated synonyms.
            organism: An Organism name or record.
            source: A `bionty.Source` record that specifies the version to validate against.
            strict_source: Determines the validation behavior against records in the registry.
                - If `False`, validation will include all records in the registry, ignoring the specified source.
                - If `True`, validation will only include records in the registry  that are linked to the specified source.
                Note: this parameter won't affect validation against public sources.

        Returns:
            If `return_mapper` is `False`: a list of standardized names. Otherwise,
            a dictionary of mapped values with mappable synonyms as keys and
            standardized names as values.

        See Also:
            :meth:`~lamindb.models.CanCurate.add_synonym`
                Add synonyms.
            :meth:`~lamindb.models.CanCurate.remove_synonym`
                Remove synonyms.

        Example::

            import bionty as bt

            # save some gene records
            bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol", organism="human").save()

            # standardize gene synonyms
            gene_synonyms = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            bt.Gene.standardize(gene_synonyms)
            #> ['A1CF', 'A1BG', 'BRCA2', 'FANCD20']
        """
        return _standardize(
            cls=cls,
            values=values,
            field=field,
            return_field=return_field,
            return_mapper=return_mapper,
            case_sensitive=case_sensitive,
            mute=mute,
            strict_source=strict_source,
            source_aware=source_aware,
            keep=keep,
            synonyms_field=synonyms_field,
            organism=organism,
            source=source,
        )

    def add_synonym(
        self,
        synonym: str | ListLike,
        force: bool = False,
        save: bool | None = None,
    ):
        """Add synonyms to a record.

        Args:
            synonym: The synonyms to add to the record.
            force: Whether to add synonyms even if they are already synonyms of other records.
            save: Whether to save the record to the database.

        See Also:
            :meth:`~lamindb.models.CanCurate.remove_synonym`
                Remove synonyms.

        Example::

            import bionty as bt

            # save "T cell" record
            record = bt.CellType.from_source(name="T cell").save()
            record.synonyms
            #> "T-cell|T lymphocyte|T-lymphocyte"

            # add a synonym
            record.add_synonym("T cells")
            record.synonyms
            #> "T cells|T-cell|T-lymphocyte|T lymphocyte"
        """
        _check_synonyms_field_exist(self)
        _add_or_remove_synonyms(
            synonym=synonym, record=self, force=force, action="add", save=save
        )

    def remove_synonym(self, synonym: str | ListLike):
        """Remove synonyms from a record.

        Args:
            synonym: The synonym values to remove.

        See Also:
            :meth:`~lamindb.models.CanCurate.add_synonym`
                Add synonyms

        Example::

            import bionty as bt

            # save "T cell" record
            record = bt.CellType.from_source(name="T cell").save()
            record.synonyms
            #> "T-cell|T lymphocyte|T-lymphocyte"

            # remove a synonym
            record.remove_synonym("T-cell")
            record.synonyms
            #> "T lymphocyte|T-lymphocyte"
        """
        _check_synonyms_field_exist(self)
        _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")

    def set_abbr(self, value: str):
        """Set value for abbr field and add to synonyms.

        Args:
            value: A value for an abbreviation.

        See Also:
            :meth:`~lamindb.models.CanCurate.add_synonym`

        Example::

            import bionty as bt

            # save an experimental factor record
            scrna = bt.ExperimentalFactor.from_source(name="single-cell RNA sequencing").save()
            assert scrna.abbr is None
            assert scrna.synonyms == "single-cell RNA-seq|single-cell transcriptome sequencing|scRNA-seq|single cell RNA sequencing"

            # set abbreviation
            scrna.set_abbr("scRNA")
            assert scrna.abbr == "scRNA"
            # synonyms are updated
            assert scrna.synonyms == "scRNA|single-cell RNA-seq|single cell RNA sequencing|single-cell transcriptome sequencing|scRNA-seq"
        """
        self.abbr = value

        if hasattr(self, "name") and value == self.name:
            pass
        else:
            try:
                self.add_synonym(value, save=False)
            except Exception as e:  # pragma: no cover
                logger.debug(
                    f"Encountered an Exception while attempting to add synonyms.\n{e}"
                )

        if not self._state.adding:  # type: ignore
            self.save()  # type: ignore
