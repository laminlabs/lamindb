from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING, Any, Literal, Union

from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from lamindb_setup import logger
from lamindb_setup.core import colors

from lamindb.base.utils import strict_classmethod

from ..errors import ValidationError
from ._from_values import (
    _format_values,
    _from_values,
    get_organism_record_from_field,
)
from .sqlrecord import SQLRecord, get_name_field

if TYPE_CHECKING:
    from collections.abc import Iterable

    import numpy as np
    import pandas as pd

    from lamindb.base.types import ListLike, StrField

    from .query_set import SQLRecordList


def map_synonyms(
    df: pd.DataFrame,
    identifiers: Iterable,
    field: str,
    *,
    case_sensitive: bool = False,
    return_mapper: bool = False,
    mute: bool = False,
    synonyms_field: str = "synonyms",
    sep: str = "|",
    keep: Literal["first", "last", False] = "first",
    mute_warning: bool = False,
) -> dict[str, str] | list[str]:
    """Maps input identifiers against a field with synonym fallback.

    Implements a three-tier matching priority:
    1. Exact case-sensitive field match (preserves original casing)
    2. Case-insensitive field match (when case_sensitive=False)
    3. Synonym match (with optional case-insensitive matching)

    Args:
        df: Reference DataFrame.
        identifiers: Identifiers that will be mapped against a field.
        field: The field representing the identifiers.
        case_sensitive: Whether the mapping is case sensitive.
        return_mapper: If True, returns {input : standardized field name}.
        mute: If True, suppresses logging of mapping statistics.
        synonyms_field: The field representing the concatenated synonyms.
        sep: Separator used to split synonyms.
        keep: {'first', 'last', False}, default 'first'
            When a synonym maps to multiple standardized values, determines
            which duplicates to mark as `pandas.DataFrame.duplicated`.
            - "first": returns the first mapped standardized value
            - "last": returns the last mapped standardized value
            - False: returns all mapped standardized values
        mute_warning: If True, suppresses warnings about list values when keep=False.

    Returns:
        - If return_mapper is False: a list of mapped field values in input order.
        - If return_mapper is True: a dictionary mapping input identifiers to
          standardized field values (only includes entries that were mapped).
    """
    import pandas as pd

    identifiers = list(identifiers)
    n_input = len(identifiers)

    # Handle empty inputs
    if (
        df.shape[0] == 0
        or n_input == 0
        or synonyms_field is None
        or synonyms_field == "None"
    ):
        return {} if return_mapper else identifiers

    # Validate inputs
    if field not in df.columns:
        raise KeyError(
            f"field '{field}' is invalid! Available fields are: {list(df.columns)}"
        )
    if synonyms_field not in df.columns:
        raise KeyError(
            f"synonyms_field '{synonyms_field}' is invalid! Available fields are: {list(df.columns)}"
        )
    if field == synonyms_field:
        raise KeyError("synonyms_field must be different from field!")

    # Track None positions before pandas converts them to NaN (pandas 3.0 + PyArrow)
    _none_positions = [i for i, v in enumerate(identifiers) if v is None]

    # Initialize mapping dataframe
    mapped_df = pd.DataFrame({"orig_ids": identifiers})
    mapped_df["__lookup__"] = to_str(
        mapped_df["orig_ids"], case_sensitive=case_sensitive
    )
    mapped_df["mapped"] = pd.NA

    # Step 1: Try exact case-sensitive match (highest priority)
    # This preserves original casing even when case_sensitive=False
    exact_field_values = set(df[field].dropna().drop_duplicates())
    exact_matches = mapped_df["orig_ids"].isin(exact_field_values)
    mapped_df.loc[exact_matches, "mapped"] = mapped_df.loc[exact_matches, "orig_ids"]

    # Step 2: For case-insensitive mode, try case-insensitive field matching
    if not case_sensitive:
        unmapped_mask = mapped_df["mapped"].isna()
        if unmapped_mask.any():
            # Build case-insensitive field map (keeps first occurrence)
            df_field = df[[field]].dropna(subset=[field])
            df_field["__lookup__"] = to_str(df_field[field], case_sensitive=False)
            df_field = df_field.drop_duplicates(subset=["__lookup__"], keep="first")
            field_map_lower = df_field.set_index("__lookup__")[field].to_dict()

            # Apply case-insensitive field map to unmapped entries
            mapped_df.loc[unmapped_mask, "mapped"] = mapped_df.loc[
                unmapped_mask, "__lookup__"
            ].map(field_map_lower)

    # Step 3: For still-unmapped terms, check synonyms
    unmapped_mask = mapped_df["mapped"].isna()
    if unmapped_mask.any():
        unmapped_terms = set(mapped_df.loc[unmapped_mask, "__lookup__"])

        syn_map = _build_synonym_map(
            df=df,
            synonyms_field=synonyms_field,
            field=field,
            unmapped_terms=unmapped_terms,
            case_sensitive=case_sensitive,
            keep=keep,
            sep=sep,
        )

        if syn_map:
            mapped_df.loc[unmapped_mask, "mapped"] = mapped_df.loc[
                unmapped_mask, "__lookup__"
            ].map(syn_map)

    # Log mapping statistics (only count actual changes, not exact matches)
    if keep is False:
        changed_mask = (~mapped_df["mapped"].isna()) & (
            mapped_df.apply(lambda row: row["mapped"] != row["orig_ids"], axis=1)
        )
    else:
        changed_mask = (~mapped_df["mapped"].isna()) & (
            mapped_df["mapped"] != mapped_df["orig_ids"]
        )
    n_mapped = changed_mask.sum()
    if n_mapped > 0 and not mute:
        s = "" if n_mapped == 1 else "s"
        logger.info(f"standardized {n_mapped}/{n_input} term{s}")

    # Return results
    if return_mapper:
        return _build_mapper(mapped_df, keep, mute_warning)
    else:
        result = _build_result_list(mapped_df, keep, mute_warning)
        # Restore None for originally-None inputs (pandas 3.0 PyArrow coerces None → NaN)
        for i in _none_positions:
            result[i] = None
        return result


def _build_synonym_map(
    df: pd.DataFrame,
    synonyms_field: str,
    field: str,
    unmapped_terms: set,
    case_sensitive: bool,
    keep: Literal["first", "last", False],
    sep: str,
) -> dict:
    """Build a synonym mapping dictionary for unmapped terms."""
    syn_series = explode_aggregated_column_to_map(
        df=df,
        agg_col=synonyms_field,
        target_col=field,
        keep=keep,
        sep=sep,
    )

    if not case_sensitive:
        # Convert synonym keys to lowercase for matching
        syn_series.index = syn_series.index.str.lower()
        # Remove duplicate synonym keys (keep first occurrence)
        syn_series = syn_series[~syn_series.index.duplicated(keep="first")]

    # Only keep synonym mappings for unmapped terms
    return {k: v for k, v in syn_series.to_dict().items() if k in unmapped_terms}


def _build_mapper(
    mapped_df: pd.DataFrame,
    keep: Literal["first", "last", False],
    mute_warning: bool,
) -> dict:
    """Build the mapper dictionary from mapped dataframe."""
    mapper_df = mapped_df[~mapped_df["mapped"].isna()].copy()
    mapper = dict(zip(mapper_df["orig_ids"], mapper_df["mapped"], strict=False))
    # Only include entries where mapping changed the value
    mapper = {k: v for k, v in mapper.items() if k != v}

    if keep is False:
        if not mute_warning:
            logger.warning(
                "returning mapper might contain lists as values when 'keep=False'"
            )
        return {
            k: v[0] if isinstance(v, list) and len(v) == 1 else v
            for k, v in mapper.items()
        }
    return mapper


def _build_result_list(
    mapped_df: pd.DataFrame,
    keep: Literal["first", "last", False],
    mute_warning: bool,
) -> list:
    """Build the result list from mapped dataframe."""
    import pandas as pd

    result = [
        m if not (m is None or (not isinstance(m, list) and pd.isna(m))) else o
        for m, o in zip(mapped_df["mapped"], mapped_df["orig_ids"], strict=False)
    ]

    if keep is False:
        if not mute_warning:
            logger.warning("returning list might contain lists when 'keep=False'")
        return [v[0] if isinstance(v, list) and len(v) == 1 else v for v in result]
    return result


def to_str(
    series_values: pd.Series | pd.Index | pd.Categorical,
    case_sensitive: bool = False,
) -> pd.Series:
    """Convert Pandas Series values to strings with case sensitive option."""
    if series_values.dtype.name == "category":
        try:
            categorical = series_values.cat
        except AttributeError:
            categorical = series_values
        if "" not in categorical.categories:
            values = categorical.add_categories("")
        else:
            values = series_values
        values = values.infer_objects().fillna("").astype(str)
    else:
        values = series_values.infer_objects().fillna("")
    if case_sensitive is False:
        values = values.str.lower()
    return values


def not_empty_none_na(values: Iterable) -> pd.Series:
    """Return values that are not empty string, None or NA."""
    import pandas as pd

    series = (
        pd.Series(values) if not isinstance(values, (pd.Series, pd.Index)) else values
    )

    return series[pd.Series(series).infer_objects().fillna("").astype(bool)]


def explode_aggregated_column_to_map(
    df,
    agg_col: str,
    target_col: str,
    keep: Literal["first", "last", False] = "first",
    sep: str = "|",
) -> pd.Series:
    """Explode values from an aggregated DataFrame column to map to a target column.

    Args:
        df: A DataFrame containing the agg_col and target_col.
        agg_col: The name of the aggregated column
        target_col: the name of the target column
        keep : {'first', 'last', False}, default 'first'
            Determines which duplicates to mark as `pandas.DataFrame.duplicated`
        sep: Splits all values of the agg_col by this separator.

    Returns:
        A pandas.Series indexed by the split values from the aggregated column
    """
    df = df[[target_col, agg_col]].drop_duplicates().dropna(subset=[agg_col])

    # subset to df with only non-empty strings in the agg_col
    df = df.loc[not_empty_none_na(df[agg_col]).index]

    df[agg_col] = df[agg_col].str.split(sep)
    df_explode = df.explode(agg_col)
    # remove rows with same values in agg_col and target_col
    df_explode = df_explode[df_explode[agg_col] != df_explode[target_col]]

    # group by the agg_col and return based on keep for the target_col values
    gb = df_explode.groupby(agg_col)[target_col]
    if keep == "first":
        return gb.first()
    elif keep == "last":
        return gb.last()
    elif keep is False:
        return gb.apply(list)
    else:
        raise ValueError(f"Invalid value for keep: {keep}")


class InspectResult:
    """Result of inspect.

    An InspectResult object of calls such as :meth:`~lamindb.models.CanCurate.inspect`.
    """

    def __init__(
        self,
        validated_df: pd.DataFrame,
        validated: list[str],
        nonvalidated: list[str],
        frac_validated: float,
        n_empty: int,
        n_unique: int,
    ) -> None:
        self._df = validated_df
        self._validated = validated
        self._non_validated = nonvalidated
        self._frac_validated = frac_validated
        self._n_empty = n_empty
        self._n_unique = n_unique
        self._synonyms_mapper: dict = {}

    @property
    def df(self) -> pd.DataFrame:
        """A DataFrame indexed by values with a boolean `__validated__` column."""
        return self._df

    @property
    def validated(self) -> list[str]:
        """List of successfully :meth:`~lamindb.models.CanCurate.validate` validated items."""
        return self._validated

    @property
    def non_validated(self) -> list[str]:
        """List of unsuccessfully :meth:`~lamindb.models.CanCurate.validate` items.

        This list can be used to remove any non-validated values such as
        genes that do not map against the specified source.
        """
        return self._non_validated

    @property
    def frac_validated(self) -> float:
        """Fraction of items that were validated."""
        return self._frac_validated

    @property
    def n_empty(self) -> int:
        """Number of empty items."""
        return self._n_empty

    @property
    def n_unique(self) -> int:
        """Number of unique items."""
        return self._n_unique

    @property
    def synonyms_mapper(self) -> dict:
        """Synonyms mapper dictionary.

        Such a dictionary maps the actual values to their synonyms
        which can be used to rename values accordingly.

        Examples:
            >>> markers = pd.DataFrame(index=["KI67","CCR7"])
            >>> synonyms_mapper = bt.CellMarker.standardize(markers.index, return_mapper=True)

            {'KI67': 'Ki67', 'CCR7': 'Ccr7'}
        """
        return self._synonyms_mapper

    def __getitem__(self, key) -> list[str]:
        """Bracket access to the inspect result."""
        if key == "validated":
            return self.validated
        elif key == "non_validated":
            return self.non_validated
        # backward compatibility below
        elif key == "mapped":
            return self.validated
        elif key == "not_mapped":
            return self.non_validated
        else:
            raise KeyError("invalid key")


def validate(
    identifiers: Iterable,
    field_values: Iterable,
    *,
    case_sensitive: bool = True,
    mute: bool = False,
    field: str | None = None,
    **kwargs,
) -> np.ndarray:
    """Check if elements in an iterable are present in a list of values.

    This function validates whether each element in `identifiers` is present in `field_values`.
    It returns a boolean numpy array indicating which elements are valid (True) or invalid (False).

    Args:
        identifiers: The iterable containing elements to be validated.
        field_values: The iterable containing valid values to check against.
        case_sensitive: If True, the comparison is case-sensitive.
        mute: If True, suppresses logging output
        field: Name of the field being validated, used in logging.
        **kwargs: Additional keyword arguments.
            logging: If provided as a boolean, overrides the 'mute' parameter.

    Returns:
        A boolean numpy array where True indicates a valid element and False an invalid one.

    Notes:
        - The function converts both `identifiers` and `field_values` to strings before comparison.
    """
    if isinstance(kwargs.get("logging"), bool):
        mute = not kwargs.get("logging")
    import pandas as pd

    _check_type_compatibility(identifiers, field_values)

    identifiers = list(identifiers)
    identifiers_idx = pd.Index(identifiers)
    identifiers_idx = to_str(identifiers_idx, case_sensitive=case_sensitive)

    field_values = to_str(field_values, case_sensitive=case_sensitive)

    # annotated what complies with the default ID
    matches = identifiers_idx.isin(field_values)
    if not mute:
        if len(identifiers) == 0:
            logger.warning("input has zero length")
        else:
            _validate_logging(
                _validate_stats(identifiers=identifiers, matches=matches), field=field
            )
    return matches


def _check_type_compatibility(identifiers: Iterable, field_values: Iterable) -> None:
    """Checks whether the identifiers and field_values have the same high level (numeric vs str/categorical) data type.

    Raises:
        TypeError: If the high level data types do not match.
    """
    import math

    import numpy as np
    import pandas as pd

    # Only look at the first element because we assume that the dtype is consistent for efficiency
    id_sample, value_sample = (
        next(iter(identifiers), None),
        next(iter(field_values), None),
    )

    def _is_nan(value) -> bool:
        if isinstance(value, (float, np.floating)):
            return math.isnan(value) or np.isnan(value)
        return False

    def _get_type_category(value):
        if isinstance(value, (int, float, complex, np.number)):
            return "numeric"
        elif isinstance(value, (str, np.str_, pd.Categorical)):
            return "str/categorical"
        return "unknown"

    # Real world data may have Nones and nan values. We can pass over them.
    if (
        id_sample is not None
        and value_sample is not None
        and not _is_nan(id_sample)
        and not _is_nan(value_sample)
    ):
        id_type, value_type = (
            _get_type_category(id_sample),
            _get_type_category(value_sample),
        )

        if id_type != value_type:
            raise TypeError(
                f"Type mismatch: identifiers are '{id_type}' but field_values are '{value_type}'."
            )


def _unique_rm_empty(idx: pd.Index):
    idx = idx.unique()
    return idx[(idx != "") & (~idx.isnull())]


def _validate_stats(identifiers: Iterable, matches: np.ndarray):
    import pandas as pd

    df_val = pd.DataFrame(data={"__validated__": matches}, index=identifiers)
    val = _unique_rm_empty(df_val.index[df_val["__validated__"]]).tolist()
    nonval = _unique_rm_empty(df_val.index[~df_val["__validated__"]]).tolist()

    n_unique = len(val) + len(nonval)
    if n_unique == 0:
        return InspectResult(
            validated_df=df_val,
            validated=val,
            nonvalidated=nonval,
            frac_validated=0,
            n_empty=0,
            n_unique=0,
        )
    n_empty = df_val.shape[0] - n_unique
    frac_nonval = round(len(nonval) / n_unique * 100, 1)
    frac_val = 100 - frac_nonval

    return InspectResult(
        validated_df=df_val,
        validated=val,
        nonvalidated=nonval,
        frac_validated=frac_val,
        n_empty=n_empty,
        n_unique=n_unique,
    )


def _validate_logging(result: InspectResult, field: str | None = None) -> None:
    """Logging of the validated result to stdout."""
    field_msg = ""
    if field is not None:
        field_msg = f" for {colors.italic(field)}"
    empty_warn_msg = ""
    if result.n_empty > 0:
        unique_s = "" if result.n_unique == 1 else "s"
        empty_s = " is" if result.n_empty == 1 else "s are"
        empty_warn_msg = (
            f"received {result.n_unique} unique term{unique_s},"
            f" {result.n_empty} empty/duplicated term{empty_s} ignored"
        )
    s = "" if len(result.validated) == 1 else "s"
    are = "is" if len(result.validated) == 1 else "are"
    success_msg = ""
    if len(result.validated) > 0:
        success_msg = (
            f"{colors.green(f'{len(result.validated)} unique term{s}')} ({result.frac_validated:.2f}%)"
            f" {are} validated{field_msg}"
        )
    if result.frac_validated < 100:
        s = "" if len(result.non_validated) == 1 else "s"
        are = "is" if len(result.non_validated) == 1 else "are"
        print_values = ", ".join([f"'{i}'" for i in result.non_validated[:10]])
        if len(result.non_validated) > 10:
            print_values += ", ..."
        warn_msg = (
            f"{colors.yellow(f'{len(result.non_validated)} unique term{s}')} ({(100 - result.frac_validated):.2f}%)"
            f" {are} not validated{field_msg}: {colors.yellow(print_values)}"
        )
        if len(empty_warn_msg) > 0:
            logger.warning(empty_warn_msg)
        if len(success_msg) > 0:
            logger.success(success_msg)
        logger.warning(warn_msg)
    else:
        logger.success(success_msg)


def inspect(
    df: pd.DataFrame,
    identifiers: Iterable,
    field: str,
    *,
    standardize: bool = True,
    mute: bool = False,
    **kwargs,
) -> InspectResult:
    """Inspect if a list of identifiers are mappable to the entity reference.

    Args:
        df: DataFrame containing the field.
        identifiers: Identifiers that will be checked against the field.
        field: The BiontyField of the ontology to compare against.
                Examples are 'ontology_id' to map against the source ID
                or 'name' to map against the ontologies field names.
        return_df: Whether to return a Pandas DataFrame.

    Returns:
        InspectResult object.
    """
    # backward compat
    if isinstance(kwargs.get("logging"), bool):
        mute = not kwargs.get("logging")
    import pandas as pd

    identifiers = list(identifiers)
    uniq_identifiers = _unique_rm_empty(pd.Index(identifiers)).tolist()
    # empty DataFrame or input
    if df.shape[0] == 0 or len(uniq_identifiers) == 0:
        result = _validate_stats(
            identifiers=identifiers,
            matches=[False] * len(identifiers),  # type:ignore
        )
        if not mute:
            _validate_logging(result=result, field=field)
        if kwargs.get("return_df") is True:
            return result.df
        else:
            return result

    # check if index is compliant with exact matches
    matches = validate(
        identifiers=identifiers, field_values=df[field], case_sensitive=True, mute=True
    )
    # matches if case sensitive is turned off
    noncs_matches = validate(
        identifiers=identifiers, field_values=df[field], case_sensitive=False, mute=True
    )

    msg_casing = "inconsistent casing/" if noncs_matches.sum() > matches.sum() else ""

    result = _validate_stats(identifiers=identifiers, matches=matches)

    # backward compat
    info_msg = ""
    if standardize and len(result.non_validated) > 0:
        try:
            synonyms_mapper = map_synonyms(
                df=df,
                identifiers=result.non_validated,
                field=field,
                return_mapper=True,
                case_sensitive=False,
                mute=True,
            )
            if len(synonyms_mapper) > 0:
                print_values = ", ".join(
                    list(synonyms_mapper.keys())[:10]  # type:ignore
                )
                if len(synonyms_mapper) > 10:
                    print_values += ", ..."
                s = "" if len(synonyms_mapper) == 1 else "s"
                labels = colors.yellow(
                    f"{len(synonyms_mapper)} unique terms with {msg_casing}synonym{s}"
                )
                info_msg = f"detected {labels}: {colors.yellow(print_values)}"
                result._synonyms_mapper = synonyms_mapper

        except Exception:  # noqa: S110
            pass
    if not mute:
        _validate_logging(result=result, field=field)
        if len(info_msg) > 0:
            logger.print(f"   {info_msg}")
            logger.print(f"→  standardize terms via {colors.italic('.standardize()')}")

    # backward compat
    if kwargs.get("return_df") is True:
        return result.df

    return result


def standardize(
    df: Any,
    identifiers: Iterable,
    field: str,
    *,
    return_field: str = None,
    case_sensitive: bool = False,
    return_mapper: bool = False,
    mute: bool = False,
    synonyms_field: str = "synonyms",
    sep: str = "|",
    keep: Literal["first", "last", False] = "first",
) -> dict[str, str] | list[str]:
    """Standardizes input identifiers against a concatenated synonyms column.

    Will also standardize casing.

    Args:
        df: Reference DataFrame.
        identifiers: Identifiers that will be mapped against a field.
        field: The field representing the identifiers.
        return_field: The field to return. Defaults to field.
        case_sensitive: Whether the mapping is case sensitive.
        return_mapper: If True, returns {input synonyms : standardized field name}.
        mute: If True, suppresses logging.
        synonyms_field: The field representing the concatenated synonyms.
        sep: Which separator is used to separate synonyms.
        keep: {'first', 'last', False}, default 'first'
            When a synonym maps to multiple standardized values, determines
            which duplicates to mark as `pandas.DataFrame.duplicated`.
            - "first": returns the first mapped standardized value
            - "last": returns the last mapped standardized value
            - False: returns all mapped standardized value

    Returns:
        - If return_mapper is False: a list of mapped field values.
        - If return_mapper is True: a dictionary of mapped values with mappable
            identifiers as keys and values mapped to field as values.
    """
    if df.shape[0] == 0 or len(identifiers) == 0:  # type: ignore
        if return_mapper:
            return {}
        else:
            return identifiers  # type: ignore

    # default return_field to field if not specified
    return_field = field if return_field is None else return_field

    # map synonyms
    result = map_synonyms(
        df=df,
        identifiers=identifiers,
        field=field,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        mute=mute,
        synonyms_field=synonyms_field,
        sep=sep,
        keep=keep,
    )

    if return_field == field:
        return result

    # convert identifiers to return_field
    # always get the full list of values (identifiers)
    if return_mapper:
        values = map_synonyms(
            df=df,
            identifiers=identifiers,
            field=field,
            return_mapper=False,
            case_sensitive=case_sensitive,
            mute=True,
            synonyms_field=synonyms_field,
            sep=sep,
            keep=keep,
            mute_warning=True,
        )
    else:
        values = result

    # no values can be converted
    if len(values) == 0:
        if not mute:
            logger.warning(
                f"no values can be converted from {field} to {return_field}!"
            )
        return values
    if keep is False:
        # flatten list of lists
        values = list(
            chain(*[item if isinstance(item, list) else [item] for item in values])
        )
    else:
        # deal with duplications here
        df = df.drop_duplicates(subset=[field], keep=keep)

    values_df = df[df[field].isin(values)]
    mapper = values_df[[field, return_field]].set_index(field)[return_field]
    if keep is False:
        mapper = (
            mapper.groupby(field)
            .agg(lambda x: list(x) if len(x) > 1 else x.iloc[0])
            .to_dict()
        )

    if return_mapper:
        # deals with the case where the mapper is a list
        return_dict: dict = {}
        for k, v in result.items():  # type: ignore
            if isinstance(v, list):
                return_dict[k] = []
                for x in v:
                    if mapper.get(x) is None:
                        continue
                    if isinstance(mapper.get(x), list):
                        return_dict[k].extend(mapper.get(x))
                    else:
                        return_dict[k].append(mapper.get(x))
            else:
                if mapper.get(v) is not None:
                    return_dict[k] = mapper.get(v)
        # add non-synonyms converted values
        return_dict.update(
            {
                k: v
                for k, v in mapper.items()
                if k
                not in set(
                    chain(*[v if isinstance(v, list) else [v] for v in result.values()])  # type: ignore
                )
            }
        )
        return return_dict
    else:
        return [mapper.get(v, v) for v in values]


def _check_if_record_in_db(record: str | SQLRecord | None, using: str | None):
    """Check if the record is from the target DB."""
    if isinstance(record, SQLRecord):
        if using is not None and using != "default":
            if record._state.db != using:
                raise ValueError(
                    f"record must be a {record.__class__.__get_name_with_module__()} record from instance '{using}'!"
                )


def _concat_lists(values: ListLike | list[list[str]] | str) -> ListLike:
    """Concatenate a list of lists of strings into a single list."""
    import pandas as pd

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
    return values  # type: ignore


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
) -> InspectResult:
    """{}"""  # noqa: D415
    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)
    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.filter().all()
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
    standardize = True
    if hasattr(registry, "_name_field") and field_str != registry._name_field:
        standardize = False

    # inspect in the DB
    result_db = inspect(
        df=_filter_queryset_with_organism(queryset=queryset, organism=organism_record),
        identifiers=values,
        field=field_str,
        standardize=standardize,
        mute=mute,
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
                standardize=standardize,
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
    import numpy as np
    import pandas as pd

    return_str = True if isinstance(values, str) else False
    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)

    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.filter().all()
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
            msg = f"Your {queryset.model.__name__} registry is empty, consider populating it first!"
            if hasattr(queryset.model, "source_id"):
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
    from_source: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    organism: str | SQLRecord | None = None,
    source: SQLRecord | None = None,
    strict_source: bool = False,
) -> list[str] | dict[str, str]:
    """{}"""  # noqa: D415
    import numpy as np
    import pandas as pd

    return_str = True if isinstance(values, str) else False
    values = _concat_lists(values)

    field_str = get_name_field(cls, field=field)
    return_field_str = get_name_field(
        cls, field=field if return_field is None else return_field
    )
    queryset = cls.all() if isinstance(cls, (QuerySet, Manager)) else cls.filter().all()
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

    # standardized names from the DB
    std_names_db = standardize(
        df=df,
        identifiers=values,
        field=field_str,
        return_field=return_field_str,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
        return_mapper=return_mapper,
        mute=mute,
    )

    def _return(result: Any, mapper: dict[str, str]):
        if return_mapper:
            return mapper
        else:
            if return_str and len(result) == 1:
                return result[0]
            return result

    # map synonyms in public source
    if hasattr(registry, "source_id") and from_source:
        mapper: dict[str, str] = {}
        if return_mapper:
            if isinstance(std_names_db, dict):
                mapper = std_names_db
            std_names_db = standardize(
                df=df,
                identifiers=values,
                field=field_str,
                return_field=return_field_str,
                case_sensitive=case_sensitive,
                keep=keep,
                synonyms_field=synonyms_field,
                return_mapper=False,
                mute=True,
            )

        val_res = registry.validate(
            std_names_db, field=field, mute=True, organism=organism_record
        )
        if all(val_res):
            return _return(result=std_names_db, mapper=mapper)

        nonval = np.array(std_names_db)[~val_res]
        std_names_bt_mapper = registry.public(
            organism=organism_record, source=source
        ).standardize(
            nonval,
            return_mapper=True,
            mute=True,
            field=field_str,
            return_field=return_field_str,
            case_sensitive=case_sensitive,
            keep=keep,
            synonyms_field=synonyms_field,
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
        # standardize() is annotated as dict | list, but a categorical Series can
        # come back when values are converted through pandas.
        names: Any = std_names_db
        if hasattr(names, "dtype") and isinstance(names.dtype, pd.CategoricalDtype):
            result = names.cat.rename_categories(std_names_bt_mapper).tolist()
        else:
            result = pd.Series(names).replace(std_names_bt_mapper).tolist()
        return _return(result=result, mapper=mapper)

    else:
        return _return(
            result=std_names_db,
            mapper=std_names_db if isinstance(std_names_db, dict) else {},
        )


def _add_or_remove_synonyms(
    synonym: str | ListLike,
    record: HasSynonyms,
    action: Literal["add", "remove"],
    force: bool = False,
    save: bool | None = None,
):
    """Add or remove synonyms."""
    assert isinstance(record, SQLRecord), "record must be a SQLRecord instance"

    def check_synonyms_in_all_records(synonyms: set[str], record: HasSynonyms):
        """Errors if input synonym is associated with other records in the DB."""
        import pandas as pd
        from IPython.display import display

        syns_all = (
            record.__class__.filter().exclude(synonyms="").exclude(synonyms=None)  # type: ignore
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

    @strict_classmethod
    def inspect(  # type: ignore[misc]
        cls: type[CanCurate],
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

        Example:

            Inspect gene symbols::

                import bionty as bt

                # populate the gene registry
                bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol", organism="human").save()

                # inspect gene symbols
                symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
                result = bt.Gene.inspect(symbols, field=bt.Gene.symbol, organism="human")
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

    @strict_classmethod
    def validate(  # type: ignore[misc]
        cls: type[CanCurate],
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

        Example:

            Validate gene symbols::

                import bionty as bt

                # populate the gene registry
                bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol", organism="human").save()

                # validate gene symbols
                symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
                bt.Gene.validate(symbols, field=bt.Gene.symbol, organism="human")
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

    @strict_classmethod
    def from_values(  # type: ignore[misc]
        cls: type[CanCurate],
        values: ListLike,
        field: StrField | None = None,
        create: bool = False,
        organism: Union[SQLRecord, str, None] = None,
        source: SQLRecord | None = None,
        mute: bool = False,
    ) -> SQLRecordList:
        """Bulk create validated records by parsing values for an identifier such as a name or an id).

        Args:
            values: A list of values for an identifier, e.g. `["name1", "name2"]`.
            field: A `SQLRecord` field to look up, e.g., `bt.CellMarker.name`.
            create: Whether to create records if they don't exist.
            organism: A `bionty.Organism` name or record.
            source: A `bionty.Source` record to validate against to create records for.
            mute: Whether to mute logging.

        Returns:
            A list of validated records. For bionty registries. Also returns knowledge-coupled records.

        Notes:
            For more info, see tutorial: :doc:`docs:manage-ontologies`.

        Example:

            Bulk create labels::

                # from invalid values logs warnings & returns an empty list
                ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"])
                assert len(ulabels) == 0

                # from valid values or via `create=True` returns label objects
                ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"], create=True).save()
                assert len(ulabels) == 3

                # bulk create cell type labels from a public ontology
                import bionty as bt
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

    @strict_classmethod
    def standardize(  # type: ignore[misc]
        cls: type[CanCurate],
        values: ListLike,
        field: StrField | None = None,
        *,
        return_field: StrField | None = None,
        return_mapper: bool = False,
        case_sensitive: bool = False,
        mute: bool = False,
        from_source: bool = True,
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
            from_source: Whether to standardize from public source. Defaults to `True` for BioRecord registries.
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
            :meth:`~lamindb.models.HasSynonyms.add_synonym`
                Add synonyms.
            :meth:`~lamindb.models.HasSynonyms.remove_synonym`
                Remove synonyms.

        Example:

            Standardize gene identifiers::

                import bionty as bt

                # save some gene objects
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
            from_source=from_source,
            keep=keep,
            synonyms_field=synonyms_field,
            organism=organism,
            source=source,
        )


class HasSynonyms:
    """Mixin for registries that define a `synonyms` field."""

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
            :meth:`~lamindb.models.HasSynonyms.remove_synonym`
                Remove synonyms.

        Example:

            Add a synonym for a cell type::

                import bionty as bt

                # create a "T cell" object
                t_cell = bt.CellType.from_source(name="T cell").save()
                t_cell.synonyms
                #> "T-cell|T lymphocyte|T-lymphocyte"

                # add a synonym
                t_cell.add_synonym("T cells")
                t_cell.synonyms
                #> "T cells|T-cell|T-lymphocyte|T lymphocyte"
        """
        _add_or_remove_synonyms(
            synonym=synonym, record=self, force=force, action="add", save=save
        )

    def remove_synonym(self, synonym: str | ListLike):
        """Remove synonyms from a record.

        Args:
            synonym: The synonym values to remove.

        See Also:
            :meth:`~lamindb.models.HasSynonyms.add_synonym`
                Add synonyms

        Example:

            Remove a synonym for a cell type::

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
        _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


class HasAbbr:
    """Mixin for registries that define an `abbr` field."""

    def set_abbr(self, value: str):
        """Set value for `abbr` field and add to `synonyms`.

        Args:
            value: A value for an abbreviation.

        See Also:
            :meth:`~lamindb.models.HasSynonyms.add_synonym`

        Example:

            Add an abbreiation for an experimental factor::

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
        elif isinstance(self, HasSynonyms):
            self.add_synonym(value, save=False)
        if not self._state.adding:  # type: ignore
            self.save()  # type: ignore
