from __future__ import annotations

import sys
from collections import defaultdict
from datetime import date, datetime  # noqa: TC003
from itertools import chain
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NamedTuple,
    overload,
)

from django.core.validators import RegexValidator
from django.db import IntegrityError, models
from django.db.models import CASCADE, PROTECT, Field, Q
from django.db.models.base import ModelBase
from django.db.models.fields.related import (
    ManyToManyField,
    ManyToManyRel,
    ManyToOneRel,
)
from lamin_utils import colors
from lamindb_setup import _check_instance_setup
from lamindb_setup.core.hashing import HASH_LENGTH, hash_dict

from lamindb.base import deprecated, doc_args
from lamindb.base.fields import (
    BigIntegerField,
    BooleanField,
    CharField,
    DateField,
    DateTimeField,
    EmailField,
    ForeignKey,
    IntegerField,
    OneToOneField,
    TextField,
    URLField,
)

from .base.ids import base62_8, base62_12, base62_20
from .base.types import (
    ArtifactKind,
    FeatureDtype,
    FieldAttr,
    ListLike,
    StrField,
    TransformType,
)
from .base.users import current_user_id

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from lamin_utils._inspect import InspectResult
    from lamindb_setup.core.types import UPathStr
    from mudata import MuData
    from pyarrow.dataset import Dataset as PyArrowDataset
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from upath import UPath

    from lamindb.core import LabelManager, MappedCollection, QuerySet, RecordList
    from lamindb.core.storage import AnnDataAccessor, BackedAccessor


_TRACKING_READY: bool | None = None


class IsVersioned(models.Model):
    """Base class for versioned models."""

    class Meta:
        abstract = True

    _len_stem_uid: int

    version: str | None = CharField(max_length=30, null=True, db_index=True)
    """Version (default `None`).

    Defines version of a family of records characterized by the same `stem_uid`.

    Consider using `semantic versioning <https://semver.org>`__
    with `Python versioning <https://peps.python.org/pep-0440/>`__.
    """
    is_latest: bool = BooleanField(default=True, db_index=True)
    """Boolean flag that indicates whether a record is the latest in its version family."""

    @overload
    def __init__(self): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        self._revises = kwargs.pop("revises") if "revises" in kwargs else None
        super().__init__(*args, **kwargs)

    @property
    def stem_uid(self) -> str:
        """Universal id characterizing the version family.

        The full uid of a record is obtained via concatenating the stem uid and version information::

            stem_uid = random_base62(n_char)  # a random base62 sequence of length 12 (transform) or 16 (artifact, collection)
            version_uid = "0000"  # an auto-incrementing 4-digit base62 number
            uid = f"{stem_uid}{version_uid}"  # concatenate the stem_uid & version_uid

        """
        return self.uid[: self._len_stem_uid]  # type: ignore

    @property
    def versions(self) -> QuerySet:
        """Lists all records of the same version family.

        >>> new_artifact = ln.Artifact(df2, revises=artifact).save()
        >>> new_artifact.versions()
        """
        db = self._state.db
        if db is not None and db != "default":
            return self.__class__.using(db).filter(uid__startswith=self.stem_uid)  # type: ignore
        else:
            return self.__class__.filter(uid__startswith=self.stem_uid)  # type: ignore

    def _add_to_version_family(self, revises: IsVersioned, version: str | None = None):
        """Add current record to a version family.

        Args:
            revises: a record that belongs to the version family.
            version: semantic version of the record.
        """
        pass


def current_run() -> Run | None:
    global _TRACKING_READY

    if not _TRACKING_READY:
        _TRACKING_READY = _check_instance_setup()
    if _TRACKING_READY:
        import lamindb.core

        return lamindb.context.run
    else:
        return None


class TracksRun(models.Model):
    """Base class tracking latest run, creating user, and `created_at` timestamp."""

    class Meta:
        abstract = True

    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""
    run: Run | None = ForeignKey(
        "lamindb.Run", PROTECT, null=True, default=current_run, related_name="+"
    )
    """Last run that created or updated the record."""

    @overload
    def __init__(self): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)


class TracksUpdates(models.Model):
    """Base class tracking previous runs and `updated_at` timestamp."""

    class Meta:
        abstract = True

    updated_at: datetime = DateTimeField(auto_now=True, db_index=True)
    """Time of last update to record."""

    @overload
    def __init__(self): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)


class CanCurate:
    """Base class providing :class:`~lamindb.core.Record`-based validation."""

    @classmethod
    def inspect(
        cls,
        values: ListLike,
        field: str | StrField | None = None,
        *,
        mute: bool = False,
        organism: str | Record | None = None,
        source: Record | None = None,
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

        See Also:
            :meth:`~lamindb.core.CanCurate.validate`

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> ln.save(bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol"))
            >>> gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            >>> result = bt.Gene.inspect(gene_symbols, field=bt.Gene.symbol)
            >>> result.validated
            ['A1CF', 'A1BG']
            >>> result.non_validated
            ['FANCD1', 'FANCD20']
        """
        pass

    @classmethod
    def validate(
        cls,
        values: ListLike,
        field: str | StrField | None = None,
        *,
        mute: bool = False,
        organism: str | Record | None = None,
        source: Record | None = None,
    ) -> np.ndarray:
        """Validate values against existing values of a string field.

        Note this is strict validation, only asserts exact matches.

        Args:
            values: Values that will be validated against the field.
            field: The field of values.
                    Examples are `'ontology_id'` to map against the source ID
                    or `'name'` to map against the ontologies field names.
            mute: Whether to mute logging.
            organism: An Organism name or record.
            source: A `bionty.Source` record that specifies the version to validate against.

        Returns:
            A vector of booleans indicating if an element is validated.

        See Also:
            :meth:`~lamindb.core.CanCurate.inspect`

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> ln.save(bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol"))
            >>> gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            >>> bt.Gene.validate(gene_symbols, field=bt.Gene.symbol)
            array([ True,  True, False, False])
        """
        pass

    def from_values(
        cls,
        values: ListLike,
        field: StrField | None = None,
        create: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
        mute: bool = False,
    ) -> RecordList:
        """Bulk create validated records by parsing values for an identifier such as a name or an id).

        Args:
            values: A list of values for an identifier, e.g.
                `["name1", "name2"]`.
            field: A `Record` field to look up, e.g., `bt.CellMarker.name`.
            create: Whether to create records if they don't exist.
            organism: A `bionty.Organism` name or record.
            source: A `bionty.Source` record to validate against to create records for.
            mute: Whether to mute logging.

        Returns:
            A list of validated records. For bionty registries. Also returns knowledge-coupled records.

        Notes:
            For more info, see tutorial: :doc:`docs:bio-registries`.

        Examples:

            Bulk create from non-validated values will log warnings & returns empty list:

            >>> ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"], field="name")
            >>> assert len(ulabels) == 0

            Bulk create records from validated values returns the corresponding existing records:

            >>> ln.save([ln.ULabel(name=name) for name in ["benchmark", "prediction", "test"]])
            >>> ulabels = ln.ULabel.from_values(["benchmark", "prediction", "test"], field="name")
            >>> assert len(ulabels) == 3

            Bulk create records from public reference:

            >>> import bionty as bt
            >>> records = bt.CellType.from_values(["T cell", "B cell"], field="name")
            >>> records
        """
        pass

    @classmethod
    def standardize(
        cls,
        values: Iterable,
        field: str | StrField | None = None,
        *,
        return_field: str | StrField | None = None,
        return_mapper: bool = False,
        case_sensitive: bool = False,
        mute: bool = False,
        public_aware: bool = True,
        keep: Literal["first", "last", False] = "first",
        synonyms_field: str = "synonyms",
        organism: str | Record | None = None,
        source: Record | None = None,
    ) -> list[str] | dict[str, str]:
        """Maps input synonyms to standardized names.

        Args:
            values: Identifiers that will be standardized.
            field: The field representing the standardized names.
            return_field: The field to return. Defaults to field.
            return_mapper: If `True`, returns `{input_value: standardized_name}`.
            case_sensitive: Whether the mapping is case sensitive.
            mute: Whether to mute logging.
            public_aware: Whether to standardize from Bionty reference. Defaults to `True` for Bionty registries.
            keep: When a synonym maps to multiple names, determines which duplicates to mark as `pd.DataFrame.duplicated`:
                    - `"first"`: returns the first mapped standardized name
                    - `"last"`: returns the last mapped standardized name
                    - `False`: returns all mapped standardized name.

                  When `keep` is `False`, the returned list of standardized names will contain nested lists in case of duplicates.

                  When a field is converted into return_field, keep marks which matches to keep when multiple return_field values map to the same field value.
            synonyms_field: A field containing the concatenated synonyms.
            organism: An Organism name or record.
            source: A `bionty.Source` record that specifies the version to validate against.

        Returns:
            If `return_mapper` is `False`: a list of standardized names. Otherwise,
            a dictionary of mapped values with mappable synonyms as keys and
            standardized names as values.

        See Also:
            :meth:`~lamindb.core.CanCurate.add_synonym`
                Add synonyms.
            :meth:`~lamindb.core.CanCurate.remove_synonym`
                Remove synonyms.

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> ln.save(bt.Gene.from_values(["A1CF", "A1BG", "BRCA2"], field="symbol"))
            >>> gene_synonyms = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
            >>> standardized_names = bt.Gene.standardize(gene_synonyms)
            >>> standardized_names
            ['A1CF', 'A1BG', 'BRCA2', 'FANCD20']
        """
        pass

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
            :meth:`~lamindb.core.CanCurate.remove_synonym`
                Remove synonyms.

        Examples:
            >>> import bionty as bt
            >>> bt.CellType.from_source(name="T cell").save()
            >>> lookup = bt.CellType.lookup()
            >>> record = lookup.t_cell
            >>> record.synonyms
            'T-cell|T lymphocyte|T-lymphocyte'
            >>> record.add_synonym("T cells")
            >>> record.synonyms
            'T cells|T-cell|T-lymphocyte|T lymphocyte'
        """
        pass

    def remove_synonym(self, synonym: str | ListLike):
        """Remove synonyms from a record.

        Args:
            synonym: The synonym values to remove.

        See Also:
            :meth:`~lamindb.core.CanCurate.add_synonym`
                Add synonyms

        Examples:
            >>> import bionty as bt
            >>> bt.CellType.from_source(name="T cell").save()
            >>> lookup = bt.CellType.lookup()
            >>> record = lookup.t_cell
            >>> record.synonyms
            'T-cell|T lymphocyte|T-lymphocyte'
            >>> record.remove_synonym("T-cell")
            'T lymphocyte|T-lymphocyte'
        """
        pass

    def set_abbr(self, value: str):
        """Set value for abbr field and add to synonyms.

        Args:
            value: A value for an abbreviation.

        See Also:
            :meth:`~lamindb.core.CanCurate.add_synonym`

        Examples:
            >>> import bionty as bt
            >>> bt.ExperimentalFactor.from_source(name="single-cell RNA sequencing").save()
            >>> scrna = bt.ExperimentalFactor.get(name="single-cell RNA sequencing")
            >>> scrna.abbr
            None
            >>> scrna.synonyms
            'single-cell RNA-seq|single-cell transcriptome sequencing|scRNA-seq|single cell RNA sequencing'
            >>> scrna.set_abbr("scRNA")
            >>> scrna.abbr
            'scRNA'
            >>> scrna.synonyms
            'scRNA|single-cell RNA-seq|single cell RNA sequencing|single-cell transcriptome sequencing|scRNA-seq'
            >>> scrna.save()
        """
        pass


class HasParents:
    """Base class for hierarchical registries (ontologies)."""

    def view_parents(
        self,
        field: StrField | None = None,
        with_children: bool = False,
        distance: int = 5,
    ):
        """View parents in an ontology.

        Args:
            field: Field to display on graph
            with_children: Whether to also show children.
            distance: Maximum distance still shown.

        Ontological hierarchies: :class:`~lamindb.ULabel` (project & sub-project), :class:`~bionty.CellType` (cell type & subtype).

        Examples:
            >>> import bionty as bt
            >>> bt.Tissue.from_source(name="subsegmental bronchus").save()
            >>> record = bt.Tissue.get(name="respiratory tube")
            >>> record.view_parents()
            >>> tissue.view_parents(with_children=True)
        """
        pass

    def query_parents(self) -> QuerySet:
        """Query parents in an ontology."""
        pass

    def query_children(self) -> QuerySet:
        """Query children in an ontology."""
        pass


class ValidateFields:
    pass


RECORD_REGISTRY_EXAMPLE = """Example::

        from lamindb import Record, fields

        # sub-classing `Record` creates a new registry
        class Experiment(Record):
            name: str = fields.CharField()

        # instantiating `Experiment` creates a record `experiment`
        experiment = Experiment(name="my experiment")

        # you can save the record to the database
        experiment.save()

        # `Experiment` refers to the registry, which you can query
        df = Experiment.filter(name__startswith="my ").df()
"""


# this is the metaclass for Record
@doc_args(RECORD_REGISTRY_EXAMPLE)
class Registry(ModelBase):
    """Metaclass for :class:`~lamindb.core.Record`.

    Each `Registry` *object* is a `Record` *class* and corresponds to a table in the metadata SQL database.

    You work with `Registry` objects whenever you use *class methods* of `Record`.

    You call any subclass of `Record` a "registry" and their objects "records". A `Record` object corresponds to a row in the SQL table.

    If you want to create a new registry, you sub-class `Record`.

    {}

    Note: `Registry` inherits from Django's `ModelBase`.
    """

    def __new__(cls, name, bases, attrs, **kwargs):
        new_class = super().__new__(cls, name, bases, attrs, **kwargs)
        return new_class

    # below creates a sensible auto-complete behavior that differs across the
    # class and instance level in Jupyter Editors it doesn't have any effect for
    # static type analyzer like pylance used in VSCode
    def __dir__(cls):
        # this is needed to bring auto-complete on the class-level back
        # https://laminlabs.slack.com/archives/C04FPE8V01W/p1717535625268849
        # Filter class attributes, excluding instance methods
        exclude_instance_methods = "sphinx" not in sys.modules
        # https://laminlabs.slack.com/archives/C04FPE8V01W/p1721134595920959

        def include_attribute(attr_name, attr_value):
            if attr_name.startswith("__"):
                return False
            if exclude_instance_methods and callable(attr_value):
                return isinstance(attr_value, (classmethod, staticmethod, type))
            return True

        # check also inherited attributes
        if hasattr(cls, "mro"):
            attrs = chain(*(c.__dict__.items() for c in cls.mro()))
        else:
            attrs = cls.__dict__.items()

        result = []
        for attr_name, attr_value in attrs:
            if attr_name not in result and include_attribute(attr_name, attr_value):
                result.append(attr_name)

        # Add non-dunder attributes from Registry
        for attr in dir(Registry):
            if not attr.startswith("__") and attr not in result:
                result.append(attr)
        return result

    def __repr__(cls) -> str:
        return registry_repr(cls)

    def lookup(
        cls,
        field: StrField | None = None,
        return_field: StrField | None = None,
    ) -> NamedTuple:
        """Return an auto-complete object for a field.

        Args:
            field: The field to look up the values for. Defaults to first string field.
            return_field: The field to return. If `None`, returns the whole record.

        Returns:
            A `NamedTuple` of lookup information of the field values with a
            dictionary converter.

        See Also:
            :meth:`~lamindb.core.Record.search`

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> bt.Gene.from_source(symbol="ADGB-DT").save()
            >>> lookup = bt.Gene.lookup()
            >>> lookup.adgb_dt
            >>> lookup_dict = lookup.dict()
            >>> lookup_dict['ADGB-DT']
            >>> lookup_by_ensembl_id = bt.Gene.lookup(field="ensembl_gene_id")
            >>> genes.ensg00000002745
            >>> lookup_return_symbols = bt.Gene.lookup(field="ensembl_gene_id", return_field="symbol")
        """
        pass

    def filter(cls, *queries, **expressions) -> QuerySet:
        """Query records.

        Args:
            queries: One or multiple `Q` objects.
            expressions: Fields and values passed as Django query expressions.

        Returns:
            A :class:`~lamindb.core.QuerySet`.

        See Also:
            - Guide: :doc:`docs:registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples:
            >>> ln.ULabel(name="my label").save()
            >>> ln.ULabel.filter(name__startswith="my").df()
        """
        pass

    def get(
        cls,
        idlike: int | str | None = None,
        **expressions,
    ) -> Record:
        """Get a single record.

        Args:
            idlike: Either a uid stub, uid or an integer id.
            expressions: Fields and values passed as Django query expressions.

        Returns:
            A record.

        Raises:
            :exc:`docs:lamindb.core.exceptions.DoesNotExist`: In case no matching record is found.

        See Also:
            - Guide: :doc:`docs:registries`
            - Django documentation: `Queries <https://docs.djangoproject.com/en/stable/topics/db/queries/>`__

        Examples:
            >>> ulabel = ln.ULabel.get("FvtpPJLJ")
            >>> ulabel = ln.ULabel.get(name="my-label")
        """
        pass

    def df(
        cls,
        include: str | list[str] | None = None,
        features: bool | list[str] = False,
        limit: int = 100,
    ) -> pd.DataFrame:
        """Convert to `pd.DataFrame`.

        By default, shows all direct fields, except `updated_at`.

        Use arguments `include` or `feature` to include other data.

        Args:
            include: Related fields to include as columns. Takes strings of
                form `"ulabels__name"`, `"cell_types__name"`, etc. or a list
                of such strings.
            features: If `True`, map all features of the
                :class:`~lamindb.Feature` registry onto the resulting
                `DataFrame`. Only available for `Artifact`.
            limit: Maximum number of rows to display from a Pandas DataFrame.
                Defaults to 100 to reduce database load.

        Examples:

            Include the name of the creator in the `DataFrame`:

            >>> ln.ULabel.df(include="created_by__name"])

            Include display of features for `Artifact`:

            >>> df = ln.Artifact.df(features=True)
            >>> ln.view(df)  # visualize with type annotations

            Only include select features:

            >>> df = ln.Artifact.df(features=["cell_type_by_expert", "cell_type_by_model"])
        """
        pass

    def search(
        cls,
        string: str,
        *,
        field: StrField | None = None,
        limit: int | None = 20,
        case_sensitive: bool = False,
    ) -> QuerySet:
        """Search.

        Args:
            string: The input string to match against the field ontology values.
            field: The field or fields to search. Search all string fields by default.
            limit: Maximum amount of top results to return.
            case_sensitive: Whether the match is case sensitive.

        Returns:
            A sorted `DataFrame` of search results with a score in column `score`.
            If `return_queryset` is `True`.  `QuerySet`.

        See Also:
            :meth:`~lamindb.core.Record.filter`
            :meth:`~lamindb.core.Record.lookup`

        Examples:
            >>> ulabels = ln.ULabel.from_values(["ULabel1", "ULabel2", "ULabel3"], field="name")
            >>> ln.save(ulabels)
            >>> ln.ULabel.search("ULabel2")
        """
        pass

    def using(
        cls,
        instance: str | None,
    ) -> QuerySet:
        """Use a non-default LaminDB instance.

        Args:
            instance: An instance identifier of form "account_handle/instance_name".

        Examples:
            >>> ln.ULabel.using("account_handle/instance_name").search("ULabel7", field="name")
                        uid    score
            name
            ULabel7  g7Hk9b2v  100.0
            ULabel5  t4Jm6s0q   75.0
            ULabel6  r2Xw8p1z   75.0
        """
        pass

    def __get_module_name__(cls) -> str:
        schema_module_name = cls.__module__.split(".")[0]
        module_name = schema_module_name.replace("lnschema_", "")
        if module_name == "lamindb":
            module_name = "core"
        return module_name

    @deprecated("__get_module_name__")
    def __get_schema_name__(cls) -> str:
        return cls.__get_module_name__()

    def __get_name_with_module__(cls) -> str:
        module_name = cls.__get_module_name__()
        if module_name == "core":
            module_prefix = ""
        else:
            module_prefix = f"{module_name}."
        return f"{module_prefix}{cls.__name__}"

    @deprecated("__get_name_with_module__")
    def __get_name_with_schema__(cls) -> str:
        return cls.__get_name_with_module__()


class BasicRecord(models.Model, metaclass=Registry):
    """Basic metadata record.

    It has the same methods as Record, but doesn't have the additional fields.

    It's mainly used for LinkORMs and similar.
    """

    class Meta:
        abstract = True


class Space(BasicRecord):
    """Spaces."""

    id: int = models.SmallAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    name: str = models.CharField(max_length=100, db_index=True)
    """Name of space."""
    description: str | None = CharField(null=True)
    """Description of space."""
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "User", CASCADE, default=None, related_name="+", null=True
    )
    """Creator of run."""


@doc_args(RECORD_REGISTRY_EXAMPLE)
class Record(BasicRecord, metaclass=Registry):
    """Metadata record.

    Every `Record` is a data model that comes with a registry in form of a SQL
    table in your database.

    Sub-classing `Record` creates a new registry while instantiating a `Record`
    creates a new record.

    {}

    `Record`'s metaclass is :class:`~lamindb.core.Registry`.

    `Record` inherits from Django's `Model` class. Why does LaminDB call it `Record`
    and not `Model`? The term `Record` can't lead to confusion with statistical,
    machine learning or biological models.
    """

    _branch_code: int = models.SmallIntegerField(db_index=True, default=1, db_default=1)
    """Whether record is on a branch, in archive or in trash.

    This dictates whether a record appears in queries & searches.

    Coding is as follows:

    - 3: template (hidden in queries & searches)
    - 2: draft (hidden in queries & searches)
    - 1: default (visible in queries & searches)
    - 0: archive (hidden, meant to be kept)
    - -1: trash (hidden, scheduled for deletion)

    Any integer higher than >3 codes a branch that's involved in a pull request.
    """
    space: Space = ForeignKey(Space, PROTECT, default=1, db_default=1)
    """The space in which the record lives."""
    aux: dict[str, Any] | None = models.JSONField(
        default=None, db_default=None, null=True
    )
    """Auxiliary field for dictionary-like metadata."""

    def save(self, *args, **kwargs) -> Record:
        """Save.

        Always saves to the default database.
        """
        # we need this here because we're using models also from plain
        # django outside of lamindb
        super().save(*args, **kwargs)
        return self

    def delete(self) -> None:
        """Delete."""
        pass

    class Meta:
        abstract = True


class FeatureManager:
    """Feature manager."""

    pass


class ParamManager:
    """Param manager."""

    pass


class ParamManagerArtifact(ParamManager):
    """Param manager."""

    pass


class ParamManagerRun(ParamManager):
    """Param manager."""

    pass


# -------------------------------------------------------------------------------------
# A note on required fields at the Record level
#
# As Django does most of its validation on the Form-level, it doesn't offer functionality
# for validating the integrity of an Record object upon instantation (similar to pydantic)
#
# For required fields, we define them as commonly done on the SQL level together
# with a validator in Record (validate_required_fields)
#
# This goes against the Django convention, but goes with the SQLModel convention
# (Optional fields can be null on the SQL level, non-optional fields cannot)
#
# Due to Django's convention where CharFieldAttr has pre-configured (null=False, default=""), marking
# a required field necessitates passing `default=None`. Without the validator it would trigger
# an error at the SQL-level, with it, it triggers it at instantiation

# -------------------------------------------------------------------------------------
# A note on class and instance methods of core Record
#
# All of these are defined and tested within lamindb, in files starting with _{orm_name}.py

# -------------------------------------------------------------------------------------
# A note on maximal lengths of char fields
#
# 100 characters:
#     "Raindrops pitter-pattered on the windowpane, blurring the"
#     "city lights outside, curled up with a mug."
# A good maximal length for a name (title).
#
# 150 characters: We choose this for name maximal length because some users like long names.
#
# 255 characters:
#     "In creating a precise 255-character paragraph, one engages in"
#     "a dance of words, where clarity meets brevity. Every syllable counts,"
#     "illustrating the skill in compact expression, ensuring the essence of the"
#     "message shines through within the exacting limit."
# This is a good maximal length for a description field.


class User(BasicRecord, CanCurate):
    """Users.

    All data in this registry is synched from `lamin.ai` to ensure a universal
    user identity. There is no need to manually create records.

    Examples:

        Query a user by handle:

        >>> user = ln.User.get(handle="testuser1")
        >>> user
    """

    _name_field: str = "handle"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=8)
    """Universal id, valid across DB instances."""
    handle: str = CharField(max_length=30, unique=True, db_index=True)
    """Universal handle, valid across DB instances (required)."""
    name: str | None = CharField(max_length=150, db_index=True, null=True)
    """Name (optional)."""  # has to match hub specification, where it's also optional
    created_artifacts: Artifact
    """Artifacts created by user."""
    created_transforms: Transform
    """Transforms created by user."""
    created_runs: Run
    """Runs created by user."""
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    updated_at: datetime = DateTimeField(auto_now=True, db_index=True)
    """Time of last update to record."""

    @overload
    def __init__(
        self,
        handle: str,
        email: str,
        name: str | None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)


class Storage(Record, TracksRun, TracksUpdates):
    """Storage locations.

    A storage location is either a directory/folder (local or in the cloud) or
    an entire S3/GCP bucket.

    A LaminDB instance can manage and link multiple storage locations. But any
    storage location is managed by *at most one* LaminDB instance.

    .. dropdown:: Managed vs. linked storage locations

        The LaminDB instance can update & delete artifacts in managed storage
        locations but merely read artifacts in linked storage locations.

        When you transfer artifacts from another instance, the default is to
        only copy metadata into the target instance, but merely link the data.

        The `instance_uid` field indicates the managing LaminDB instance of a
        storage location.

        When you delete a LaminDB instance, you'll be warned about data in managed
        storage locations while data in linked storage locations is ignored.

    See Also:
        :attr:`~lamindb.core.Settings.storage`
            Default storage.
        :attr:`~lamindb.setup.core.StorageSettings`
            Storage settings.

    Examples:

        Configure the default storage location upon initiation of a LaminDB instance::

            lamin init --storage ./mydata # or "s3://my-bucket" or "gs://my-bucket"

        View the default storage location:

        >>> ln.settings.storage
        PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata')

        Dynamically change the default storage:

        >>> ln.settings.storage = "./storage_2" # or a cloud bucket
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "root"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, max_length=12, default=base62_12, db_index=True)
    """Universal id, valid across DB instances."""
    # we are very conservative here with 255 characters
    root: str = CharField(max_length=255, db_index=True, unique=True)
    """Root path of storage. n s3 path.  local path, etc. (required)."""
    description: str | None = CharField(max_length=255, db_index=True, null=True)
    """A description of what the storage location is used for (optional)."""
    type: str = CharField(max_length=30, db_index=True)
    """Can be "local" vs. "s3" vs. "gs"."""
    region: str | None = CharField(max_length=64, db_index=True, null=True)
    """Cloud storage region, if applicable."""
    instance_uid: str | None = CharField(max_length=12, db_index=True, null=True)
    """Instance that manages this storage location."""
    artifacts: Artifact
    """Artifacts contained in this storage location."""

    @overload
    def __init__(
        self,
        root: str,
        type: str,
        region: str | None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

    @property
    def path(self) -> Path | UPath:
        """Bucket or folder path.

        Cloud storage bucket:

        >>> ln.Storage("s3://my-bucket").save()

        Directory/folder in cloud storage:

        >>> ln.Storage("s3://my-bucket/my-directory").save()

        Local directory/folder:

        >>> ln.Storage("./my-directory").save()
        """
        pass


# does not inherit from TracksRun because the Transform
# is needed to define a run
class Transform(Record, IsVersioned):
    """Data transformations.

    A "transform" can refer to a Python function, a script, a notebook, or a
    pipeline. If you execute a transform, you generate a run
    (:class:`~lamindb.Run`). A run has inputs and outputs.

    A pipeline is typically created with a workflow tool (Nextflow, Snakemake,
    Prefect, Flyte, MetaFlow, redun, Airflow, ...) and stored in a versioned
    repository.

    Transforms are versioned so that a given transform version maps on a given
    source code version.

    .. dropdown:: Can I sync transforms to git?

        If you switch on
        :attr:`~lamindb.core.Settings.sync_git_repo` a script-like transform is
        synched to its hashed state in a git repository upon calling `ln.track()`.

        >>> ln.settings.sync_git_repo = "https://github.com/laminlabs/lamindb"
        >>> ln.track()

    The definition of transforms and runs is consistent the OpenLineage
    specification where a :class:`~lamindb.Transform` record would be called a
    "job" and a :class:`~lamindb.Run` record a "run".

    Args:
        name: `str` A name or title.
        key: `str | None = None` A short name or path-like semantic key.
        type: `TransformType | None = "pipeline"` See :class:`~lamindb.base.types.TransformType`.
        revises: `Transform | None = None` An old version of the transform.

    See Also:
        :meth:`~lamindb.core.Context.track`
            Globally track a script, notebook or pipeline run.
        :class:`~lamindb.Run`
            Executions of transforms.

    Notes:
        - :doc:`docs:track`
        - :doc:`docs:data-flow`
        - :doc:`docs:redun`
        - :doc:`docs:nextflow`
        - :doc:`docs:snakemake`

    Examples:

        Create a transform for a pipeline:

        >>> transform = ln.Transform(name="Cell Ranger", version="7.2.0", type="pipeline").save()

        Create a transform from a notebook:

        >>> ln.track()

        View predecessors of a transform:

        >>> transform.view_lineage()
    """

    class Meta(Record.Meta, IsVersioned.Meta):
        abstract = False

    _len_stem_uid: int = 12
    _len_full_uid: int = 16
    _name_field: str = "key"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=_len_full_uid)
    """Universal id."""
    key: str | None = CharField(db_index=True, null=True)
    """A name or "/"-separated path-like string.

    All transforms with the same key are part of the same version family.
    """
    description: str | None = CharField(db_index=True, null=True)
    """A description."""
    type: TransformType = CharField(
        max_length=20,
        db_index=True,
        default="pipeline",
    )
    """:class:`~lamindb.base.types.TransformType` (default `"pipeline"`)."""
    source_code: str | None = TextField(null=True)
    """Source code of the transform.

    .. versionchanged:: 0.75
       The `source_code` field is no longer an artifact, but a text field.
    """
    hash: str | None = CharField(max_length=HASH_LENGTH, db_index=True, null=True)
    """Hash of the source code."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """Reference for the transform, e.g., a URL."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Reference type of the transform, e.g., 'url'."""
    runs: Run
    """Runs of this transform."""
    ulabels: ULabel = models.ManyToManyField(
        "ULabel", through="TransformULabel", related_name="transforms"
    )
    """ULabel annotations of this transform."""
    predecessors: Transform = models.ManyToManyField(
        "self", symmetrical=False, related_name="successors"
    )
    """Preceding transforms.

    These are auto-populated whenever an artifact or collection serves as a run
    input, e.g., `artifact.run` and `artifact.transform` get populated & saved.

    The table provides a more convenient method to query for the predecessors that
    bypasses querying the :class:`~lamindb.Run`.

    It also allows to manually add predecessors whose outputs are not tracked in a run.
    """
    successors: Transform
    """Subsequent transforms.

    See :attr:`~lamindb.Transform.predecessors`.
    """
    output_artifacts: Artifact
    """The artifacts generated by all runs of this transform.

    If you're looking for the outputs of a single run, see :attr:`lamindb.Run.output_artifacts`.
    """
    output_collections: Collection
    """The collections generated by all runs of this transform.

    If you're looking for the outputs of a single run, see :attr:`lamindb.Run.output_collections`.
    """
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    updated_at: datetime = DateTimeField(auto_now=True, db_index=True)
    """Time of last update to record."""
    created_by: User = ForeignKey(
        User, PROTECT, default=current_user_id, related_name="created_transforms"
    )
    """Creator of record."""
    _template: Transform | None = ForeignKey(
        "Transform", PROTECT, related_name="_derived_from", default=None, null=True
    )
    """Creating template."""

    @overload
    def __init__(
        self,
        name: str,
        key: str | None = None,
        type: TransformType | None = None,
        revises: Transform | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

    @property
    def name(self) -> str:
        """Name of the transform.

        Splits `key` on `/` and returns the last element.
        """
        return self.key.split("/")[-1]

    @property
    def latest_run(self) -> Run:
        """The latest run of this transform."""
        pass

    def view_lineage(self) -> None:
        """View lineage of transforms."""
        pass


class Param(Record, CanCurate, TracksRun, TracksUpdates):
    """Parameters of runs & models."""

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    name: str = CharField(max_length=100, db_index=True)
    dtype: str = CharField(max_length=64, db_index=True)
    """Data type ("num", "cat", "int", "float", "bool", "datetime").

    For categorical types, can define from which registry values are
    sampled, e.g., `cat[ULabel]` or `cat[bionty.CellType]`.
    """
    type: Param | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of param (e.g., 'Pipeline', 'ModelTraining', 'PostProcessing').

    Allows to group features by type, e.g., all read outs, all metrics, etc.
    """
    records: Param
    """Records of this type."""
    is_type: bool = BooleanField(default=None, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    _expect_many: bool = models.BooleanField(default=False, db_default=False)
    """Indicates whether values for this param are expected to occur a single or multiple times for an artifact/run (default `False`).

    - if it's `False` (default), the values mean artifact/run-level values and a dtype of `datetime` means `datetime`
    - if it's `True`, the values are from an aggregation, which this seems like an edge case but when characterizing a model ensemble trained with different parameters it could be relevant
    """
    # backward fields
    values: ParamValue
    """Values for this parameter."""


# FeatureValue behaves in many ways like a link in a LinkORM
# in particular, we don't want a _public field on it
# Also, we don't inherit from TracksRun because a ParamValue
# is typically created before a run is created and we want to
# avoid delete cycles (for Model params though it might be helpful)
class ParamValue(Record):
    """Parameter values.

    Is largely analogous to `FeatureValue`.
    """

    # we do not have a unique constraint on param & value because it leads to hashing errors
    # for large dictionaries: https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0000
    # we do not hash values because we have `get_or_create` logic all over the place
    # and also for checking whether the (param, value) combination exists
    # there does not seem an issue with querying for a dict-like value
    # https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0001
    _name_field: str = "value"

    param: Param = ForeignKey(Param, CASCADE, related_name="values")
    """The dimension metadata."""
    value: Any = (
        models.JSONField()
    )  # stores float, integer, boolean, datetime or dictionaries
    """The JSON-like value."""
    # it'd be confusing and hard to populate a run here because these
    # values are typically created upon creating a run
    # hence, ParamValue does _not_ inherit from TracksRun but manually
    # adds created_at & created_by
    # because ParamValue cannot be updated, we don't need updated_at
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    created_by: User = ForeignKey(
        User, PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""
    hash: str = CharField(max_length=HASH_LENGTH, null=True, db_index=True)

    class Meta:
        constraints = [
            # For simple types, use direct value comparison
            models.UniqueConstraint(
                fields=["param", "value"],
                name="unique_simple_param_value",
                condition=Q(hash__isnull=True),
            ),
            # For complex types (dictionaries), use hash
            models.UniqueConstraint(
                fields=["param", "hash"],
                name="unique_complex_param_value",
                condition=Q(hash__isnull=False),
            ),
        ]

    @classmethod
    def get_or_create(cls, param, value):
        # Simple types: int, float, str, bool
        if isinstance(value, (int, float, str, bool)):
            try:
                return cls.objects.create(param=param, value=value, hash=None), False
            except IntegrityError:
                return cls.objects.get(param=param, value=value), True

        # Complex types: dict, list
        else:
            hash = hash_dict(value)
            try:
                return cls.objects.create(param=param, value=value, hash=hash), False
            except IntegrityError:
                return cls.objects.get(param=param, hash=hash), True


class Run(Record):
    """Runs of transforms.

    Args:
        transform: `Transform` A :class:`~lamindb.Transform` record.
        reference: `str | None = None` For instance, an external ID or a download URL.
        reference_type: `str | None = None` For instance, `redun_id`, `nextflow_id` or `url`.

    See Also:
        :meth:`~lamindb.core.Context.track`
            Track global run & transform records for a notebook or pipeline.

    Examples:

        Create a run record:

        >>> ln.Transform(name="Cell Ranger", version="7.2.0", type="pipeline").save()
        >>> transform = ln.Transform.get(name="Cell Ranger", version="7.2.0")
        >>> run = ln.Run(transform)

        Create a global run context for a custom transform:

        >>> ln.track(transform=transform)
        >>> ln.context.run  # globally available run

        Track a global run context for a notebook or script:

        >>> ln.track()  # Jupyter notebook metadata is automatically parsed
        >>> ln.context.run
    """

    _name_field: str = "started_at"

    params: ParamManager = ParamManagerRun  # type: ignore
    """Param manager.

    Guide: :ref:`track-run-parameters`

    Example::

        run.params.add_values({
            "learning_rate": 0.01,
            "input_dir": "s3://my-bucket/mydataset",
            "downsample": True,
            "preprocess_params": {
                "normalization_type": "cool",
                "subset_highlyvariable": True,
            },
        })
    """

    id: int = models.BigAutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=20, default=base62_20)
    """Universal id, valid across DB instances."""
    name: str | None = CharField(max_length=150, null=True)
    """A name."""
    transform = ForeignKey(Transform, CASCADE, related_name="runs")
    """The transform :class:`~lamindb.Transform` that is being run."""
    started_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Start time of run."""
    finished_at: datetime | None = DateTimeField(db_index=True, null=True, default=None)
    """Finished time of run."""
    # we don't want to make below a OneToOne because there could be the same trivial report
    # generated for many different runs
    report: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_report_of", default=None
    )
    """Report of run, e.g.. n html file."""
    _logfile: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_logfile_of", default=None
    )
    """Report of run, e.g.. n html file."""
    environment: Artifact | None = ForeignKey(
        "Artifact", PROTECT, null=True, related_name="_environment_of", default=None
    )
    """Computational environment for the run.

    For instance, `Dockerfile`, `docker image`, `requirements.txt`, `environment.yml`, etc.
    """
    input_artifacts: Artifact
    """The artifacts serving as input for this run.

    Related accessor: :attr:`~lamindb.Artifact.input_of_runs`.
    """
    output_artifacts: Artifact
    """The artifacts generated by this run.

    Related accessor: via :attr:`~lamindb.Artifact.run`
    """
    input_collections: Collection
    """The collections serving as input for this run."""
    output_collections: Collection
    """The collections generated by this run."""
    _param_values: ParamValue = models.ManyToManyField(
        ParamValue, through="RunParamValue", related_name="runs"
    )
    """Parameter values."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like a URL or external ID (such as from a workflow manager)."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of reference such as a workflow manager execution ID."""
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of first creation. Mismatches ``started_at`` if the run is re-run."""
    created_by: User = ForeignKey(
        User, CASCADE, default=current_user_id, related_name="created_runs"
    )
    """Creator of run."""
    initiated_by_run: Run | None = ForeignKey(
        "Run", CASCADE, null=True, related_name="initiated_runs", default=None
    )
    """The run that triggered the current run.

    This is not a preceding run. The preceding runs ("predecessors") is the set
    of runs that produced the output artifacts that serve as the inputs for the
    present run.

    Be careful with using this field at this point.
    """
    children: Run
    """The runs that are triggered by this run."""
    _is_consecutive: bool | None = BooleanField(null=True)
    """Indicates whether code was consecutively executed. Is relevant for notebooks."""
    _status_code: int = models.SmallIntegerField(default=0, db_index=True)
    """Status code of the run.

    - 0: scheduled
    - 1: started
    - 2: errored
    - 3: aborted
    - 4: completed
    """

    @overload
    def __init__(
        self,
        transform: Transform,
        reference: str | None = None,
        reference_type: str | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)


class ULabel(Record, HasParents, CanCurate, TracksRun, TracksUpdates):
    """Universal labels.

    Args:
        name: `str` A name.
        description: `str` A description.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    A `ULabel` record provides the easiest way to annotate a dataset
    with a label: `"My project"`, `"curated"`, or `"Batch X"`:

        >>> my_project = ULabel(name="My project")
        >>> my_project.save()
        >>> artifact.ulabels.add(my_project)

    Often, a ulabel is measured *within* a dataset. For instance, an artifact
    might characterize 2 species of the Iris flower (`"setosa"` &
    `"versicolor"`) measured by a `"species"` feature. Use the
    :class:`~lamindb.Curator` flow to automatically parse, validate, and
    annotate with labels that are contained in `DataFrame` or `AnnData`
    artifacts.

    .. note::

        If you work with complex entities like cell lines, cell types, tissues,
        etc., consider using the pre-defined biological registries in
        :mod:`bionty` to label artifacts & collections.

        If you work with biological samples, likely, the only sustainable way of
        tracking metadata, is to create a custom schema module.

    See Also:
        :meth:`~lamindb.Feature`
            Dimensions of measurement for artifacts & collections.
        :attr:`~lamindb.Artifact.features`
            Feature manager for an artifact.

    Examples:

        Create a new label:

        >>> train_split = ln.ULabel(name="train").save()

        Organize labels in a hierarchy:

        >>> split_type = ln.ULabel(name="Split", is_type=True).save()
        >>> train_split = ln.ULabel(name="train", type="split_type").save()

        Label an artifact:

        >>> artifact.ulabels.add(ulabel)

        Query by `ULabel`:

        >>> ln.Artifact.filter(ulabels=train_split)
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=8, default=base62_8)
    """A universal random id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True)
    """Name or title of ulabel (`unique=True`)."""
    type: ULabel | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of ulabel, e.g., `"donor"`, `"split"`, etc.

    Allows to group ulabels by type, e.g., all donors, all split ulabels, etc.
    """
    records: ULabel
    """Records of this type."""
    is_type: bool = BooleanField(default=None, db_index=True, null=True)
    """Distinguish types from instances of the type.

    For example, a ulabel "Project" would be a type, and the actual projects "Project 1", "Project 2", would be records of that `type`.
    """
    description: str | None = TextField(null=True)
    """A description (optional)."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like URL or external ID."""
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of reference such as a donor_id from Vendor X."""
    parents: ULabel = models.ManyToManyField(
        "self", symmetrical=False, related_name="children"
    )
    """Parent entities of this ulabel.

    For advanced use cases, you can build an ontology under a given `type`.

    Say, if you modeled `CellType` as a `ULabel`, you would introduce a type `CellType` and model the hiearchy of cell types under it.
    """
    children: ULabel
    """Child entities of this ulabel.

    Reverse accessor for parents.
    """
    transforms: Transform
    """Transforms annotated with this ulabel."""
    artifacts: Artifact
    """Artifacts annotated with this ulabel."""
    collections: Collection
    """Collections annotated with this ulabel."""

    @overload
    def __init__(
        self,
        name: str,
        description: str | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        pass


class Feature(Record, CanCurate, TracksRun, TracksUpdates):
    """Dataset dimensions.

    A feature represents a dimension of a dataset, such as a column in a
    `DataFrame`. The `Feature` registry organizes metadata of features.

    The `Feature` registry helps you organize and query datasets based on their
    features and corresponding label annotations. For instance, when working
    with a "T cell" label, it could be measured through different features
    such as `"cell_type_by_expert"` where an expert manually classified the
    cell, or `"cell_type_by_model"` where a computational model made the
    classification.

    The two most important metadata of a feature are its `name` and the `dtype`.
    In addition to typical data types, LaminDB has a `"num"` `dtype` to
    concisely denote the union of all numerical types.

    Args:
        name: `str` Name of the feature, typically.  column name.
        dtype: `FeatureDtype | Registry | list[Registry]` See :class:`~lamindb.base.types.FeatureDtype`.
            For categorical types, can define from which registry values are
            sampled, e.g., `ULabel` or `[ULabel, bionty.CellType]`.
        unit: `str | None = None` Unit of measure, ideally SI (`"m"`, `"s"`, `"kg"`, etc.) or `"normalized"` etc.
        description: `str | None = None` A description.
        synonyms: `str | None = None` Bar-separated synonyms.

    Note:

        For more control, you can use :mod:`bionty` registries to manage simple
        biological entities like genes, proteins & cell markers. Or you define
        custom registries to manage high-level derived features like gene sets.

    See Also:
        :meth:`~lamindb.Feature.from_df`
            Create feature records from DataFrame.
        :attr:`~lamindb.Artifact.features`
            Feature manager of an artifact or collection.
        :class:`~lamindb.ULabel`
            Universal labels.
        :class:`~lamindb.FeatureSet`
            Feature sets.

    Example:

        A simple `"str"` feature.

        >>> ln.Feature(
        ...     name="sample_note",
        ...     dtype="str",
        ... ).save()

        A dtype `"cat[ULabel]"` can be more easily passed as below.

        >>> ln.Feature(
        ...     name="project",
        ...     dtype=ln.ULabel,
        ... ).save()

        A dtype `"cat[ULabel|bionty.CellType]"` can be more easily passed as below.

        >>> ln.Feature(
        ...     name="cell_type",
        ...     dtype=[ln.ULabel, bt.CellType],
        ... ).save()

    Hint:

        *Features* and *labels* denote two ways of using entities to organize data:

        1. A feature qualifies *what* is measured, i.e., a numerical or categorical random variable
        2. A label *is* a measured value, i.e., a category

        Consider annotating a dataset by that it measured expression of 30k
        genes: genes relate to the dataset as feature identifiers through a
        feature set with 30k members. Now consider annotating the artifact by
        whether that it measured the knock-out of 3 genes: here, the 3 genes act
        as labels of the dataset.

        Re-shaping data can introduce ambiguity among features & labels. If this
        happened, ask yourself what the joint measurement was: a feature
        qualifies variables in a joint measurement. The canonical data matrix
        lists jointly measured variables in the columns.

    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=12, default=base62_12)
    """Universal id, valid across DB instances."""
    name: str = CharField(max_length=150, db_index=True, unique=True)
    """Name of feature (`unique=True`)."""
    dtype: FeatureDtype = CharField(max_length=64, db_index=True)
    """Data type (:class:`~lamindb.base.types.FeatureDtype`).

    For categorical types, can define from which registry values are
    sampled, e.g., `'cat[ULabel]'` or `'cat[bionty.CellType]'`. Unions are also
    allowed if the feature samples from two registries, e.g., `'cat[ULabel|bionty.CellType]'`
    """
    type: Feature | None = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of feature (e.g., 'Readout', 'Metric', 'Metadata', 'ExpertAnnotation', 'ModelPrediction').

    Allows to group features by type, e.g., all read outs, all metrics, etc.
    """
    records: Feature
    """Records of this type."""
    is_type: bool = BooleanField(default=None, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    unit: str | None = CharField(max_length=30, db_index=True, null=True)
    """Unit of measure, ideally SI (`m`, `s`, `kg`, etc.) or 'normalized' etc. (optional)."""
    description: str | None = TextField(db_index=True, null=True)
    """A description."""
    synonyms: str | None = TextField(null=True)
    """Bar-separated (|) synonyms (optional)."""
    # we define the below ManyToMany on the feature model because it parallels
    # how other registries (like Gene, Protein, etc.) relate to FeatureSet
    # it makes the API more consistent
    feature_sets: FeatureSet = models.ManyToManyField(
        "FeatureSet", through="FeatureSetFeature", related_name="features"
    )
    """Feature sets linked to this feature."""
    _expect_many: bool = models.BooleanField(default=True, db_default=True)
    """Indicates whether values for this feature are expected to occur a single or multiple times for an artifact (default `True`).

    - if it's `True` (default), the values come from an observation-level aggregation and a dtype of `datetime` on the observation-level mean `set[datetime]` on the artifact-level
    - if it's `False` it's an artifact-level value and datetime means datetime; this is an edge case because an arbitrary artifact would always be a set of arbitrary measurements that would need to be aggregated ("one just happens to measure a single cell line in that artifact")
    """
    # backward fields
    values: FeatureValue
    """Values for this feature."""

    @overload
    def __init__(
        self,
        name: str,
        dtype: FeatureDtype | Registry | list[Registry],
        unit: str | None,
        description: str | None,
        synonyms: str | None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        pass

    @classmethod
    def from_df(cls, df: pd.DataFrame, field: FieldAttr | None = None) -> RecordList:
        """Create Feature records for columns."""
        pass

    def save(self, *args, **kwargs) -> Feature:
        """Save."""
        pass


class FeatureValue(Record, TracksRun):
    """Non-categorical features values.

    Categorical feature values are stored in their respective registries:
    :class:`~lamindb.ULabel`, :class:`~bionty.CellType`, etc.

    Unlike for ULabel, in `FeatureValue`, values are grouped by features and
    not by an ontological hierarchy.
    """

    # we do not have a unique constraint on feature & value because it leads to hashing errors
    # for large dictionaries: https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0000
    # we do not hash values because we have `get_or_create` logic all over the place
    # and also for checking whether the (feature, value) combination exists
    # there does not seem an issue with querying for a dict-like value
    # https://lamin.ai/laminlabs/lamindata/transform/jgTrkoeuxAfs0001

    _name_field: str = "value"

    feature: Feature | None = ForeignKey(
        Feature, CASCADE, null=True, related_name="values", default=None
    )
    """The dimension metadata."""
    value: Any = models.JSONField()
    """The JSON-like value."""
    hash: str = CharField(max_length=HASH_LENGTH, null=True, db_index=True)
    """Value hash."""

    class Meta(BasicRecord.Meta, TracksRun.Meta):
        constraints = [
            # For simple types, use direct value comparison
            models.UniqueConstraint(
                fields=["feature", "value"],
                name="unique_simple_feature_value",
                condition=Q(hash__isnull=True),
            ),
            # For complex types (dictionaries), use hash
            models.UniqueConstraint(
                fields=["feature", "hash"],
                name="unique_complex_feature_value",
                condition=Q(hash__isnull=False),
            ),
        ]

    @classmethod
    def get_or_create(cls, feature, value):
        # Simple types: int, float, str, bool
        if isinstance(value, (int, float, str, bool)):
            try:
                return cls.objects.create(
                    feature=feature, value=value, hash=None
                ), False
            except IntegrityError:
                return cls.objects.get(feature=feature, value=value), True

        # Complex types: dict, list
        else:
            hash = hash_dict(value)
            try:
                return cls.objects.create(
                    feature=feature, value=value, hash=hash
                ), False
            except IntegrityError:
                return cls.objects.get(feature=feature, hash=hash), True


class FeatureSet(Record, CanCurate, TracksRun):
    """Feature sets.

    Stores references to sets of :class:`~lamindb.Feature` and other registries
    that may be used to identify features (e.g., :class:`~bionty.Gene` or
    :class:`~bionty.Protein`).

    .. dropdown:: Why does LaminDB model feature sets, not just features?

        1. Performance: Imagine you measure the same panel of 20k transcripts in
           1M samples. By modeling the panel as a feature set, you can link all
           your artifacts against one feature set and only need to store 1M
           instead of 1M x 20k = 20B links.
        2. Interpretation: Model protein panels, gene panels, etc.
        3. Data integration: Feature sets provide the information that determines whether two datasets can be meaningfully concatenated.

        These reasons do not hold for label sets. Hence, LaminDB does not model label sets.

    Args:
        features: `Iterable[Record]` An iterable of :class:`~lamindb.Feature`
            records to hash, e.g., `[Feature(...), Feature(...)]`. Is turned into
            a set upon instantiation. If you'd like to pass values, use
            :meth:`~lamindb.FeatureSet.from_values` or
            :meth:`~lamindb.FeatureSet.from_df`.
        dtype: `str | None = None` The simple type. Defaults to
            `None` for sets of :class:`~lamindb.Feature` records.
            Otherwise defaults to `"num"` (e.g., for sets of :class:`~bionty.Gene`).
        name: `str | None = None` A name.

    Note:

        A feature set can be identified by the `hash` its feature uids.
        It's stored in the `.hash` field.

        A `slot` provides a string key to access feature sets.
        It's typically the accessor within the registered data object, here `pd.DataFrame.columns`.


    See Also:
        :meth:`~lamindb.FeatureSet.from_values`
            Create from values.
        :meth:`~lamindb.FeatureSet.from_df`
            Create from dataframe columns.

    Examples:

        Create a featureset from df with types:

        >>> df = pd.DataFrame({"feat1": [1, 2], "feat2": [3.1, 4.2], "feat3": ["cond1", "cond2"]})
        >>> feature_set = ln.FeatureSet.from_df(df)

        Create a featureset from features:

        >>> features = [ln.Feature(name=feat, dtype="float").save() for feat in ["feat1", "feat2"]]
        >>> feature_set = ln.FeatureSet(features)

        Create a featureset from feature values:

        >>> import bionty as bt
        >>> feature_set = ln.FeatureSet.from_values(adata.var["ensemble_id"], Gene.ensembl_gene_id, organism="mouse")
        >>> feature_set.save()

        Link a feature set to an artifact:

        >>> artifact.features.add_feature_set(feature_set, slot="var")

        Link features to an artifact (will create a featureset under the hood):

        >>> artifact.features.add_values(features)
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _name_field: str = "name"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=20)
    """A universal id (hash of the set of feature values)."""
    name: str | None = CharField(max_length=150, null=True)
    """A name."""
    n = IntegerField()
    """Number of features in the set."""
    dtype: str | None = CharField(max_length=64, null=True)
    """Data type, e.g., "num", "float", "int". Is `None` for :class:`~lamindb.Feature`.

    For :class:`~lamindb.Feature`, types are expected to be heterogeneous and defined on a per-feature level.
    """
    type: Feature | None = ForeignKey(
        "self", PROTECT, null=True, related_name="records"
    )
    """Type of feature set (e.g., 'ExpressionPanel', 'ProteinPanel', 'Multimodal', 'Metadata', 'Embedding').

    Allows to group feature sets by type, e.g., all meassurements evaluating gene expression vs. protein expression vs. multi modal.
    """
    records: Feature
    """Records of this type."""
    is_type: bool = BooleanField(default=None, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    registry: str = CharField(max_length=120, db_index=True)
    """The registry that stores the feature identifiers, e.g., `'Feature'` or `'bionty.Gene'`.

    Depending on the registry, `.members` stores, e.g. `Feature` or `Gene` records.
    """
    hash: str | None = CharField(
        max_length=HASH_LENGTH, db_index=True, null=True, unique=True
    )
    """The hash of the set."""
    features: Feature
    """The features related to a `FeatureSet` record."""
    artifacts: Artifact
    """The artifacts related to a `FeatureSet` record."""

    @overload
    def __init__(
        self,
        features: Iterable[Record],
        dtype: str | None = None,
        name: str | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        pass

    @classmethod
    def from_values(  # type: ignore
        cls,
        values: ListLike,
        field: FieldAttr = Feature.name,
        type: str | None = None,
        name: str | None = None,
        mute: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
        raise_validation_error: bool = True,
    ) -> FeatureSet:
        """Create feature set for validated features.

        Args:
            values: A list of values, like feature names or ids.
            field: The field of a reference registry to map values.
            type: The simple type.
                Defaults to `None` if reference registry is :class:`~lamindb.Feature`,
                defaults to `"float"` otherwise.
            name: A name.
            organism: An organism to resolve gene mapping.
            source: A public ontology to resolve feature identifier mapping.
            raise_validation_error: Whether to raise a validation error if some values are not valid.

        Raises:
            ValidationError: If some values are not valid.

        Examples:

            >>> features = [ln.Feature(name=feat, dtype="str").save() for feat in ["feat11", "feat21"]]
            >>> feature_set = ln.FeatureSet.from_values(features)

            >>> genes = ["ENSG00000139618", "ENSG00000198786"]
            >>> feature_set = ln.FeatureSet.from_values(features, bt.Gene.ensembl_gene_id, "float")
        """
        pass

    @classmethod
    def from_df(
        cls,
        df: pd.DataFrame,
        field: FieldAttr = Feature.name,
        name: str | None = None,
        mute: bool = False,
        organism: Record | str | None = None,
        source: Record | None = None,
    ) -> FeatureSet | None:
        """Create feature set for validated features."""
        pass

    def save(self, *args, **kwargs) -> FeatureSet:
        """Save."""
        pass

    @property
    def members(self) -> QuerySet:
        """A queryset for the individual records of the set."""
        pass


class Artifact(Record, IsVersioned, TracksRun, TracksUpdates):
    """Datasets & models stored as files, folders, or arrays.

    Artifacts manage data in local or remote storage.

    Some artifacts are array-like, e.g., when stored as `.parquet`, `.h5ad`,
    `.zarr`, or `.tiledb`.

    Args:
        data: `UPathStr` A path to a local or remote folder or file.
        type: `Literal["dataset", "model"] | None = None` The artifact type.
        key: `str | None = None` A path-like key to reference artifact in default storage, e.g., `"myfolder/myfile.fcs"`. Artifacts with the same key form a revision family.
        description: `str | None = None` A description.
        revises: `Artifact | None = None` Previous version of the artifact. Triggers a revision.
        run: `Run | None = None` The run that creates the artifact.

    .. dropdown:: Typical storage formats & their API accessors

        Arrays:

        - Table: `.csv`, `.tsv`, `.parquet`, `.ipc` ⟷ `DataFrame`, `pyarrow.Table`
        - Annotated matrix: `.h5ad`, `.h5mu`, `.zrad` ⟷ `AnnData`, `MuData`
        - Generic array: HDF5 group, zarr group, TileDB store ⟷ HDF5, zarr, TileDB loaders

        Non-arrays:

        - Image: `.jpg`, `.png` ⟷ `np.ndarray`, ...
        - Fastq: `.fastq` ⟷ /
        - VCF: `.vcf` ⟷ /
        - QC: `.html` ⟷ /

        You'll find these values in the `suffix` & `accessor` fields.

        LaminDB makes some default choices (e.g., serialize a `DataFrame` as a `.parquet` file).

    See Also:
        :class:`~lamindb.Storage`
            Storage locations for artifacts.
        :class:`~lamindb.Collection`
            Collections of artifacts.
        :meth:`~lamindb.Artifact.from_df`
            Create an artifact from a `DataFrame`.
        :meth:`~lamindb.Artifact.from_anndata`
            Create an artifact from an `AnnData`.

    Examples:

        Create an artifact from a file path and pass `description`:

        >>> artifact = ln.Artifact("s3://my_bucket/my_folder/my_file.csv", description="My file")
        >>> artifact = ln.Artifact("./my_local_file.jpg", description="My image")

        You can also pass `key` to create a virtual filepath hierarchy:

        >>> artifact = ln.Artifact("./my_local_file.jpg", key="example_datasets/dataset1.jpg")

        What works for files also works for folders:

        >>> artifact = ln.Artifact("s3://my_bucket/my_folder", description="My folder")
        >>> artifact = ln.Artifact("./my_local_folder", description="My local folder")
        >>> artifact = ln.Artifact("./my_local_folder", key="project1/my_target_folder")

        .. dropdown:: Why does the API look this way?

            It's inspired by APIs building on AWS S3.

            Both boto3 and quilt select a bucket (akin to default storage in LaminDB) and define a target path through a `key` argument.

            In `boto3 <https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3/bucket/upload_file.html>`__::

                # signature: S3.Bucket.upload_file(filepath, key)
                import boto3
                s3 = boto3.resource('s3')
                bucket = s3.Bucket('mybucket')
                bucket.upload_file('/tmp/hello.txt', 'hello.txt')

            In `quilt3 <https://docs.quiltdata.com/api-reference/bucket>`__::

                # signature: quilt3.Bucket.put_file(key, filepath)
                import quilt3
                bucket = quilt3.Bucket('mybucket')
                bucket.put_file('hello.txt', '/tmp/hello.txt')


        Make a new version of an artifact:

        >>> artifact = ln.Artifact.from_df(df, key="example_datasets/dataset1.parquet").save()
        >>> artifact_v2 = ln.Artifact(df_updated, key="example_datasets/dataset1.parquet").save()

        Alternatively, if you don't want to provide a value for `key`, you can use `revises`:

        >>> artifact = ln.Artifact.from_df(df, description="My dataframe").save()
        >>> artifact_v2 = ln.Artifact(df_updated, revises=artifact).save()

    """

    class Meta(Record.Meta, IsVersioned.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _len_full_uid: int = 20
    _len_stem_uid: int = 16

    params: ParamManager = ParamManagerArtifact  # type: ignore
    """Param manager.

    Example::

        artifact.params.add_values({
            "hidden_size": 32,
            "bottleneck_size": 16,
            "batch_size": 32,
            "preprocess_params": {
                "normalization_type": "cool",
                "subset_highlyvariable": True,
            },
        })
    """

    features: FeatureManager = FeatureManager  # type: ignore
    """Feature manager.

    Features denote dataset dimensions, i.e., the variables that measure labels & numbers.

    Annotate with features & values::

       artifact.features.add_values({
            "species": organism,  # here, organism is an Organism record
            "scientist": ['Barbara McClintock', 'Edgar Anderson'],
            "temperature": 27.6,
            "study": "Candidate marker study"
       })

    Query for features & values::

        ln.Artifact.features.filter(scientist="Barbara McClintock")

    Features may or may not be part of the artifact content in storage. For
    instance, the :class:`~lamindb.Curator` flow validates the columns of a
    `DataFrame`-like artifact and annotates it with features corresponding to
    these columns. `artifact.features.add_values`, by contrast, does not
    validate the content of the artifact.
    """

    @property
    def labels(self) -> LabelManager:
        """Label manager.

        To annotate with labels, you typically use the registry-specific accessors,
        for instance :attr:`~lamindb.Artifact.ulabels`::

            candidate_marker_study = ln.ULabel(name="Candidate marker study").save()
            artifact.ulabels.add(candidate_marker_study)

        Similarly, you query based on these accessors::

            ln.Artifact.filter(ulabels__name="Candidate marker study").all()

        Unlike the registry-specific accessors, the `.labels` accessor provides
        a way of associating labels with features::

            study = ln.Feature(name="study", dtype="cat").save()
            artifact.labels.add(candidate_marker_study, feature=study)

        Note that the above is equivalent to::

            artifact.features.add_values({"study": candidate_marker_study})
        """
        from lamindb.core._label_manager import LabelManager

        return LabelManager(self)

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, db_index=True, max_length=_len_full_uid)
    """A universal random id."""
    key: str | None = CharField(db_index=True, null=True)
    """A (virtual) relative file path within the artifact's storage location.

    Setting a `key` is useful to automatically group artifacts into a version family.

    LaminDB defaults to a virtual file path to make renaming of data in object storage easy.

    If you register existing files in a storage location, the `key` equals the
    actual filepath on the underyling filesytem or object store.
    """
    description: str | None = CharField(db_index=True, null=True)
    """A description.

    LaminDB doesn't require you to pass a key, you can
    """
    storage: Storage = ForeignKey(Storage, PROTECT, related_name="artifacts")
    """Storage location, e.g. an S3 or GCP bucket or a local directory."""
    suffix: str = CharField(max_length=30, db_index=True)
    # Initially, we thought about having this be nullable to indicate folders
    # But, for instance, .zarr is stored in a folder that ends with a .zarr suffix
    """Path suffix or empty string if no canonical suffix exists.

    This is either a file suffix (`".csv"`, `".h5ad"`, etc.) or the empty string "".
    """
    kind: ArtifactKind | None = CharField(
        max_length=20,
        db_index=True,
        null=True,
    )
    """:class:`~lamindb.base.types.ArtifactKind` (default `None`)."""
    otype: str | None = CharField(max_length=64, db_index=True, null=True)
    """Default Python object type, e.g., DataFrame, AnnData."""
    size: int | None = BigIntegerField(null=True, db_index=True, default=None)
    """Size in bytes.

    Examples: 1KB is 1e3 bytes, 1MB is 1e6, 1GB is 1e9, 1TB is 1e12 etc.
    """
    hash: str | None = CharField(max_length=HASH_LENGTH, db_index=True, null=True)
    """Hash or pseudo-hash of artifact content.

    Useful to ascertain integrity and avoid duplication.
    """
    n_files: int | None = BigIntegerField(null=True, db_index=True, default=None)
    """Number of files for folder-like artifacts, `None` for file-like artifacts.

    Note that some arrays are also stored as folders, e.g., `.zarr` or `.tiledbsoma`.

    .. versionchanged:: 1.0
        Renamed from `n_objects` to `n_files`.
    """
    n_observations: int | None = BigIntegerField(null=True, db_index=True, default=None)
    """Number of observations.

    Typically, this denotes the first array dimension.
    """
    _hash_type: str | None = CharField(max_length=30, db_index=True, null=True)
    """Type of hash."""
    ulabels: ULabel = models.ManyToManyField(
        ULabel, through="ArtifactULabel", related_name="artifacts"
    )
    """The ulabels measured in the artifact (:class:`~lamindb.ULabel`)."""
    run: Run | None = ForeignKey(
        Run, PROTECT, related_name="output_artifacts", null=True, default=None
    )
    """Run that created the artifact."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_artifacts")
    """Runs that use this artifact as an input."""
    # if the artifact is replicated or update in a new run, we link the previous
    # run in previous_runs
    _previous_runs: Run = models.ManyToManyField(
        "Run", related_name="_output_artifacts_with_later_updates"
    )
    """Sequence of runs that created or updated the record."""
    collections: Collection
    """The collections that this artifact is part of."""
    feature_sets: FeatureSet = models.ManyToManyField(
        FeatureSet, related_name="artifacts", through="ArtifactFeatureSet"
    )
    """The feature sets measured in the artifact."""
    _feature_values: FeatureValue = models.ManyToManyField(
        FeatureValue, through="ArtifactFeatureValue", related_name="artifacts"
    )
    """Non-categorical feature values for annotation."""
    _param_values: ParamValue = models.ManyToManyField(
        ParamValue, through="ArtifactParamValue", related_name="artifacts"
    )
    """Parameter values."""
    _key_is_virtual: bool = BooleanField()
    """Indicates whether `key` is virtual or part of an actual file path."""
    # be mindful that below, passing related_name="+" leads to errors
    _actions: Artifact = models.ManyToManyField(
        "self", symmetrical=False, related_name="_action_targets"
    )
    """Actions to attach for the UI."""
    created_by: User = ForeignKey(
        "lamindb.User",
        PROTECT,
        default=current_user_id,
        related_name="created_artifacts",
    )
    """Creator of record."""
    _curator: dict[str, str] | None = models.JSONField(
        default=None, db_default=None, null=True
    )
    _overwrite_versions: bool = BooleanField(default=None)
    """Indicates whether to store or overwrite versions.

    It defaults to False for file-like artifacts and to True for folder-like artifacts.
    """

    @overload
    def __init__(
        self,
        # we're not choosing the name "path" for this arg because
        # it'd be confusing with `artifact.path`, which is not the same
        # so "data" conveys better that this is input data that's ingested
        # and will be moved to a target path at `artifact.path`
        # also internally, we sometimes pass "data objects" like a DataFrame
        # here; and we might refactor this but we might also keep that internal
        # usage
        data: UPathStr,
        type: ArtifactKind | None = None,
        key: str | None = None,
        description: str | None = None,
        revises: Artifact | None = None,
        run: Run | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        pass

    @property
    def transform(self) -> Transform | None:
        """Transform whose run created the artifact."""
        return self.run.transform if self.run is not None else None

    @property
    @deprecated("n_files")
    def n_objects(self) -> int:
        return self.n_files

    @property
    def path(self) -> Path:
        """Path.

        File in cloud storage, here AWS S3:

        >>> artifact = ln.Artifact("s3://my-bucket/my-file.csv").save()
        >>> artifact.path
        S3Path('s3://my-bucket/my-file.csv')

        File in local storage:

        >>> ln.Artifact("./myfile.csv", description="myfile").save()
        >>> artifact = ln.Artifact.get(description="myfile")
        >>> artifact.path
        PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata/myfile.csv')
        """
        pass

    @classmethod
    def from_df(
        cls,
        df: pd.DataFrame,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from `DataFrame`, validate & link features.

        Args:
            df: A `DataFrame` object.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.parquet"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.

        See Also:
            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Examples:
            >>> df = ln.core.datasets.df_iris_in_meter_batch1()
            >>> df.head()
              sepal_length sepal_width petal_length petal_width iris_organism_code
            0        0.051       0.035        0.014       0.002                 0
            1        0.049       0.030        0.014       0.002                 0
            2        0.047       0.032        0.013       0.002                 0
            3        0.046       0.031        0.015       0.002                 0
            4        0.050       0.036        0.014       0.002                 0
            >>> artifact = ln.Artifact.from_df(df, description="Iris flower collection batch1")
            >>> artifact.save()
        """
        pass

    @classmethod
    def from_anndata(
        cls,
        adata: AnnData | UPathStr,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from ``AnnData``, validate & link features.

        Args:
            adata: An `AnnData` object or a path of AnnData-like.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.h5ad"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.

        See Also:

            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> adata = ln.core.datasets.anndata_with_obs()
            >>> artifact = ln.Artifact.from_anndata(adata, description="mini anndata with obs")
            >>> artifact.save()
        """
        pass

    @classmethod
    def from_mudata(
        cls,
        mdata: MuData,
        key: str | None = None,
        description: str | None = None,
        run: Run | None = None,
        revises: Artifact | None = None,
        **kwargs,
    ) -> Artifact:
        """Create from ``MuData``, validate & link features.

        Args:
            mdata: An `MuData` object.
            key: A relative path within default storage,
                e.g., `"myfolder/myfile.h5mu"`.
            description: A description.
            revises: An old version of the artifact.
            run: The run that creates the artifact.

        See Also:
            :meth:`~lamindb.Collection`
                Track collections.
            :class:`~lamindb.Feature`
                Track features.

        Examples:
            >>> import bionty as bt
            >>> bt.settings.organism = "human"
            >>> mdata = ln.core.datasets.mudata_papalexi21_subset()
            >>> artifact = ln.Artifact.from_mudata(mdata, description="a mudata object")
            >>> artifact.save()
        """
        pass

    @classmethod
    def from_dir(
        cls,
        path: UPathStr,
        key: str | None = None,
        *,
        run: Run | None = None,
    ) -> list[Artifact]:
        """Create a list of artifact objects from a directory.

        Hint:
            If you have a high number of files (several 100k) and don't want to
            track them individually, create a single :class:`~lamindb.Artifact` via
            ``Artifact(path)`` for them. See, e.g., :doc:`docs:rxrx`.

        Args:
            path: Source path of folder.
            key: Key for storage destination. If `None` and
                directory is in a registered location, the inferred `key` will
                reflect the relative position. If `None` and directory is outside
                of a registered storage location, the inferred key defaults to `path.name`.
            run: A `Run` object.

        Examples:
            >>> dir_path = ln.core.datasets.generate_cell_ranger_files("sample_001", ln.settings.storage)
            >>> artifacts = ln.Artifact.from_dir(dir_path)
            >>> ln.save(artifacts)
        """
        pass

    def replace(
        self,
        data: UPathStr,
        run: Run | None = None,
        format: str | None = None,
    ) -> None:
        """Replace artifact content.

        Args:
            data: A file path.
            run: The run that created the artifact gets
                auto-linked if ``ln.track()`` was called.

        Examples:
            Say we made a change to the content of an artifact, e.g., edited the image
            `paradisi05_laminopathic_nuclei.jpg`.

            This is how we replace the old file in storage with the new file:

            >>> artifact.replace("paradisi05_laminopathic_nuclei.jpg")
            >>> artifact.save()

            Note that this neither changes the storage key nor the filename.

            However, it will update the suffix if it changes.
        """
        pass

    def open(
        self, mode: str = "r", is_run_input: bool | None = None
    ) -> (
        AnnDataAccessor
        | BackedAccessor
        | SOMACollection
        | SOMAExperiment
        | PyArrowDataset
    ):
        """Return a cloud-backed data object.

        Works for `AnnData` (`.h5ad` and `.zarr`), generic `hdf5` and `zarr`,
        `tiledbsoma` objects (`.tiledbsoma`), `pyarrow` compatible formats.

        Args:
            mode: can only be `"w"` (write mode) for `tiledbsoma` stores,
                otherwise should be always `"r"` (read-only mode).

        Notes:
            For more info, see tutorial: :doc:`/arrays`.

        Examples:

            Read AnnData in backed mode from cloud:

            >>> artifact = ln.Artifact.get(key="lndb-storage/pbmc68k.h5ad")
            >>> artifact.open()
            AnnDataAccessor object with n_obs × n_vars = 70 × 765
                constructed for the AnnData object pbmc68k.h5ad
                ...
        """
        pass

    def load(self, is_run_input: bool | None = None, **kwargs) -> Any:
        """Cache and load into memory.

        See all :mod:`~lamindb.core.loaders`.

        Examples:

            Load a `DataFrame`-like artifact:

            >>> artifact.load().head()
            sepal_length sepal_width petal_length petal_width iris_organism_code
            0        0.051       0.035        0.014       0.002                 0
            1        0.049       0.030        0.014       0.002                 0
            2        0.047       0.032        0.013       0.002                 0
            3        0.046       0.031        0.015       0.002                 0
            4        0.050       0.036        0.014       0.002                 0

            Load an `AnnData`-like artifact:

            >>> artifact.load()
            AnnData object with n_obs × n_vars = 70 × 765

            Fall back to :meth:`~lamindb.Artifact.cache` if no in-memory representation is configured:

            >>> artifact.load()
            PosixPath('/home/runner/work/lamindb/lamindb/docs/guide/mydata/.lamindb/jb7BY5UJoQVGMUOKiLcn.jpg')
        """
        pass

    def cache(self, is_run_input: bool | None = None) -> Path:
        """Download cloud artifact to local cache.

        Follows synching logic: only caches an artifact if it's outdated in the local cache.

        Returns a path to a locally cached on-disk object (say a `.jpg` file).

        Examples:

            Sync file from cloud and return the local path of the cache:

            >>> artifact.cache()
            PosixPath('/home/runner/work/Caches/lamindb/lamindb-ci/lndb-storage/pbmc68k.h5ad')
        """
        pass

    def delete(
        self, permanent: bool | None = None, storage: bool | None = None
    ) -> None:
        """Trash or permanently delete.

        A first call to `.delete()` puts an artifact into the trash (sets `_branch_code` to `-1`).
        A second call permanently deletes the artifact.

        FAQ: :doc:`docs:faq/storage`

        Args:
            permanent: Permanently delete the artifact (skip trash).
            storage: Indicate whether you want to delete the artifact in storage.

        Examples:

            For an `Artifact` object `artifact`, call:

            >>> artifact.delete()
        """
        pass

    def save(self, upload: bool | None = None, **kwargs) -> Artifact:
        """Save to database & storage.

        Args:
            upload: Trigger upload to cloud storage in instances with hybrid storage mode.

        Examples:
            >>> artifact = ln.Artifact("./myfile.csv", description="myfile")
            >>> artifact.save()
        """
        pass

    def restore(self) -> None:
        """Restore from trash.

        Examples:

            For any `Artifact` object `artifact`, call:

            >>> artifact.restore()
        """
        pass

    def describe(self) -> None:
        """Describe relations of record.

        Examples:
            >>> artifact.describe()
        """
        pass


class Collection(Record, IsVersioned, TracksRun, TracksUpdates):
    """Collections of artifacts.

    Collections provide a simple way of versioning collections of artifacts.

    Args:
        artifacts: `list[Artifact]` A list of artifacts.
        name: `str` A name.
        description: `str | None = None` A description.
        revises: `Collection | None = None` An old version of the collection.
        run: `Run | None = None` The run that creates the collection.
        meta: `Artifact | None = None` An artifact that defines metadata for the collection.
        reference: `str | None = None` For instance, an external ID or a URL.
        reference_type: `str | None = None` For instance, `"url"`.

    See Also:
        :class:`~lamindb.Artifact`

    Examples:

        Create a collection from a list of :class:`~lamindb.Artifact` objects:

        >>> collection = ln.Collection([artifact1, artifact2], name="My collection")

        Create a collection that groups a data & a metadata artifact (e.g., here :doc:`docs:rxrx`):

        >>> collection = ln.Collection(data_artifact, name="My collection", meta=metadata_artifact)

    """

    class Meta(Record.Meta, IsVersioned.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    _len_full_uid: int = 20
    _len_stem_uid: int = 16
    _name_field: str = "key"

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(
        unique=True, db_index=True, max_length=_len_full_uid, default=base62_20
    )
    """Universal id, valid across DB instances."""
    key: str = CharField(db_index=True)
    """Name or path-like key."""
    description: str | None = TextField(null=True)
    """A description or title."""
    hash: str | None = CharField(max_length=HASH_LENGTH, db_index=True, null=True)
    """Hash of collection content. 86 base64 chars allow to store 64 bytes, 512 bits."""
    reference: str | None = CharField(max_length=255, db_index=True, null=True)
    """A reference like URL or external ID."""
    # also for reference_type here, we allow an extra long max_length
    reference_type: str | None = CharField(max_length=25, db_index=True, null=True)
    """Type of reference, e.g., cellxgene Census collection_id."""
    ulabels: ULabel = models.ManyToManyField(
        "ULabel", through="CollectionULabel", related_name="collections"
    )
    """ULabels sampled in the collection (see :class:`~lamindb.Feature`)."""
    run: Run | None = ForeignKey(
        Run, PROTECT, related_name="output_collections", null=True, default=None
    )
    """:class:`~lamindb.Run` that created the `collection`."""
    input_of_runs: Run = models.ManyToManyField(Run, related_name="input_collections")
    """Runs that use this collection as an input."""
    _previous_runs: Run = models.ManyToManyField(
        "Run", related_name="_output_collections_with_later_updates"
    )
    """Sequence of runs that created or updated the record."""
    artifacts: Artifact = models.ManyToManyField(
        "Artifact", related_name="collections", through="CollectionArtifact"
    )
    """Artifacts in collection."""
    meta_artifact: Artifact | None = OneToOneField(
        "Artifact",
        PROTECT,
        null=True,
        unique=True,
        related_name="_meta_of_collection",
    )
    """An artifact that stores metadata that indexes a collection.

    It has a 1:1 correspondence with an artifact. If needed, you can access the
    collection from the artifact via a private field:
    `artifact._meta_of_collection`.
    """
    _actions: Artifact = models.ManyToManyField(Artifact, related_name="+")
    """Actions to attach for the UI."""

    @overload
    def __init__(
        self,
        artifacts: list[Artifact],
        name: str,
        description: str | None = None,
        meta: Any | None = None,
        reference: str | None = None,
        reference_type: str | None = None,
        run: Run | None = None,
        revises: Collection | None = None,
    ): ...

    @overload
    def __init__(
        self,
        *db_args,
    ): ...

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        pass

    def append(self, artifact: Artifact, run: Run | None = None) -> Collection:
        """Add an artifact to the collection.

        Creates a new version of the collection.

        Args:
            artifact: An artifact to add to the collection.
            run: The run that creates the new version of the collection.

        .. versionadded:: 0.76.14
        """
        pass

    def mapped(
        self,
        layers_keys: str | list[str] | None = None,
        obs_keys: str | list[str] | None = None,
        obsm_keys: str | list[str] | None = None,
        obs_filter: tuple[str, str | tuple[str, ...]] | None = None,
        join: Literal["inner", "outer"] | None = "inner",
        encode_labels: bool | list[str] = True,
        unknown_label: str | dict[str, str] | None = None,
        cache_categories: bool = True,
        parallel: bool = False,
        dtype: str | None = None,
        stream: bool = False,
        is_run_input: bool | None = None,
    ) -> MappedCollection:
        """Return a map-style dataset.

        Returns a `pytorch map-style dataset
        <https://pytorch.org/docs/stable/data.html#map-style-datasets>`__ by
        virtually concatenating `AnnData` arrays.

        If your `AnnData` collection is in the cloud, move them into a local
        cache first via :meth:`~lamindb.Collection.cache`.

        `__getitem__` of the `MappedCollection` object takes a single integer index
        and returns a dictionary with the observation data sample for this index from
        the `AnnData` objects in the collection. The dictionary has keys for `layers_keys`
        (`.X` is in `"X"`), `obs_keys`, `obsm_keys` (under `f"obsm_{key}"`) and also `"_store_idx"`
        for the index of the `AnnData` object containing this observation sample.

        .. note::

            For a guide, see :doc:`docs:scrna-mappedcollection`.

            This method currently only works for collections of `AnnData` artifacts.

        Args:
            layers_keys: Keys from the ``.layers`` slot. ``layers_keys=None`` or ``"X"`` in the list
                retrieves ``.X``.
            obs_keys: Keys from the ``.obs`` slots.
            obsm_keys: Keys from the ``.obsm`` slots.
            obs_filter: Select only observations with these values for the given obs column.
                Should be a tuple with an obs column name as the first element
                and filtering values (a string or a tuple of strings) as the second element.
            join: `"inner"` or `"outer"` virtual joins. If ``None`` is passed,
                does not join.
            encode_labels: Encode labels into integers.
                Can be a list with elements from ``obs_keys``.
            unknown_label: Encode this label to -1.
                Can be a dictionary with keys from ``obs_keys`` if ``encode_labels=True``
                or from ``encode_labels`` if it is a list.
            cache_categories: Enable caching categories of ``obs_keys`` for faster access.
            parallel: Enable sampling with multiple processes.
            dtype: Convert numpy arrays from ``.X``, ``.layers`` and ``.obsm``
            stream: Whether to stream data from the array backend.
            is_run_input: Whether to track this collection as run input.

        Examples:
            >>> import lamindb as ln
            >>> from torch.utils.data import DataLoader
            >>> ds = ln.Collection.get(description="my collection")
            >>> mapped = collection.mapped(obs_keys=["cell_type", "batch"])
            >>> dl = DataLoader(mapped, batch_size=128, shuffle=True)
        """
        pass

    def cache(self, is_run_input: bool | None = None) -> list[UPath]:
        """Download cloud artifacts in collection to local cache.

        Follows synching logic: only caches outdated artifacts.

        Returns paths to locally cached on-disk artifacts.

        Args:
            is_run_input: Whether to track this collection as run input.
        """
        pass

    def load(
        self,
        join: Literal["inner", "outer"] = "outer",
        is_run_input: bool | None = None,
        **kwargs,
    ) -> Any:
        """Stage and load to memory.

        Returns in-memory representation if possible such as a concatenated `DataFrame` or `AnnData` object.
        """
        pass

    def delete(self, permanent: bool | None = None) -> None:
        """Delete collection.

        Args:
            permanent: Whether to permanently delete the collection record (skips trash).

        Examples:

            For any `Collection` object `collection`, call:

            >>> collection.delete()
        """
        pass

    def save(self, using: str | None = None) -> Collection:
        """Save the collection and underlying artifacts to database & storage.

        Args:
            using: The database to which you want to save.

        Examples:
            >>> collection = ln.Collection("./myfile.csv", name="myfile")
            >>> collection.save()
        """
        pass

    def restore(self) -> None:
        """Restore collection record from trash.

        Examples:

            For any `Collection` object `collection`, call:

            >>> collection.restore()
        """
        pass

    @property
    def transform(self) -> Transform | None:
        """Transform whose run created the collection."""
        return self.run.transform if self.run is not None else None

    @property
    def name(self) -> str:
        """Name of the collection.

        Splits `key` on `/` and returns the last element.
        """
        return self.key.split("/")[-1]

    @property
    def ordered_artifacts(self) -> QuerySet:
        """Ordered `QuerySet` of `.artifacts`.

        Accessing the many-to-many field `collection.artifacts` directly gives
        you non-deterministic order.

        Using the property `.ordered_artifacts` allows to iterate through a set
        that's ordered in the order of creation.
        """
        pass

    @property
    def data_artifact(self) -> Artifact | None:
        """Access to a single data artifact.

        If the collection has a single data & metadata artifact, this allows access via::

           collection.data_artifact  # first & only element of collection.artifacts
           collection.meta_artifact  # metadata

        """
        pass

    def describe(self) -> None:
        """Describe relations of record.

        Examples:
            >>> artifact.describe()
        """
        pass


# -------------------------------------------------------------------------------------
# Project management


class Person(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """Persons.

    This registry is distinct from `User` and purely exists for project management.

    You'll soon be able to conveniently create persons from users.

    Example:
        >>> person = Person(
        ...     name="Jane Doe",
        ...     email="jane.doe@example.com",
        ...     internal=True,
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, max_length=8, db_index=True, default=base62_8)
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Name of the person (forename(s) lastname)."""
    email: str | None = EmailField(null=True, default=None)
    """Email of the person."""
    external: bool = BooleanField(default=True, db_index=True)
    """Whether the person is external to the organization."""


class Project(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """Projects.

    Example:
        >>> project = Project(
        ...     name="My Project Name",
        ...     abbr="MPN",
        ...     url="https://example.com/my_project",
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, max_length=12, db_index=True, default=base62_12)
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the Project."""
    abbr: str | None = CharField(max_length=32, db_index=True, unique=True, null=True)
    """A unique abbreviation."""
    url: str | None = URLField(max_length=255, null=True, default=None)
    """A URL to view."""
    contributors: Person = models.ManyToManyField(Person, related_name="projects")
    """Contributors associated with this project."""
    references: Reference = models.ManyToManyField("Reference", related_name="projects")
    """References associated with this project."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactProject", related_name="projects"
    )
    """Artifacts labeled with this Project."""
    transforms: Transform = models.ManyToManyField(
        Transform, through="TransformProject", related_name="projects"
    )
    """Transforms labeled with this project."""
    collections: Collection = models.ManyToManyField(
        Collection, through="CollectionProject", related_name="projects"
    )
    """Collections labeled with this project."""


class Reference(Record, CanCurate, TracksRun, TracksUpdates, ValidateFields):
    """References such as internal studies, papers, documents, or URLs.

    Example:
        >>> reference = Reference(
        ...     name="A Paper Title",
        ...     abbr="APT",
        ...     url="https://doi.org/10.1000/xyz123",
        ...     pubmed_id=12345678,
        ...     doi="10.1000/xyz123",
        ...     preprint=False,
        ...     journal="Nature Biotechnology",
        ...     description="A groundbreaking research paper.",
        ...     text="A really informative abstract.",
        ...     published_at=date(2023, 11, 21),
        ... ).save()
    """

    class Meta(Record.Meta, TracksRun.Meta, TracksUpdates.Meta):
        abstract = False

    id: int = models.AutoField(primary_key=True)
    """Internal id, valid only in one DB instance."""
    uid: str = CharField(unique=True, max_length=12, db_index=True, default=base62_12)
    """Universal id, valid across DB instances."""
    name: str = CharField(db_index=True)
    """Title or name of the reference document."""
    abbr: str | None = CharField(
        max_length=32,
        db_index=True,
        unique=True,
        null=True,
    )
    """A unique abbreviation for the reference."""
    type: Param | None = ForeignKey("self", PROTECT, null=True, related_name="records")
    """Type of reference (e.g., 'Study', 'Paper', 'Preprint').

    Allows to group reference by type, e.g., internal studies vs. all papers etc.
    """
    records: Param
    """Records of this type."""
    is_type: bool = BooleanField(default=None, db_index=True, null=True)
    """Distinguish types from instances of the type."""
    url: str | None = URLField(null=True)
    """URL linking to the reference."""
    pubmed_id: int | None = BigIntegerField(null=True, db_index=True)
    """A PudMmed ID."""
    doi: str | None = CharField(
        null=True,
        db_index=True,
        validators=[
            RegexValidator(
                regex=r"^(?:https?://(?:dx\.)?doi\.org/|doi:|DOI:)?10\.\d+/.*$",
                message="Must be a DOI (e.g., 10.1000/xyz123 or https://doi.org/10.1000/xyz123)",
            )
        ],
    )
    """Digital Object Identifier (DOI) for the reference."""
    preprint: bool = BooleanField(default=False, db_index=True)
    """Whether the reference is from a preprint."""
    public: bool = BooleanField(default=True, db_index=True)
    """Whether the reference is public."""
    journal: str | None = TextField(null=True)
    """Name of the journal."""
    description: str | None = TextField(null=True)
    """Description of the reference."""
    text: str | None = TextField(null=True)
    """Abstract or full text of the reference."""
    published_at: date | None = DateField(null=True, default=None)
    """Publication date."""
    authors: Person = models.ManyToManyField(Person, related_name="references")
    """All people associated with this reference."""
    artifacts: Artifact = models.ManyToManyField(
        Artifact, through="ArtifactReference", related_name="references"
    )
    """Artifacts labeled with this reference."""
    transforms: Artifact = models.ManyToManyField(
        Transform, through="TransformReference", related_name="references"
    )
    """Transforms labeled with this reference."""
    collections: Artifact = models.ManyToManyField(
        Collection, through="CollectionReference", related_name="references"
    )
    """Collections labeled with this reference."""


# -------------------------------------------------------------------------------------
# Link models


class LinkORM:
    pass


class FeatureSetFeature(BasicRecord, LinkORM):
    id: int = models.BigAutoField(primary_key=True)
    # we follow the lower() case convention rather than snake case for link models
    featureset: FeatureSet = ForeignKey(FeatureSet, CASCADE, related_name="+")
    feature: Feature = ForeignKey(Feature, PROTECT, related_name="+")

    class Meta:
        unique_together = ("featureset", "feature")


class ArtifactFeatureSet(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_feature_set")
    # we follow the lower() case convention rather than snake case for link models
    featureset: FeatureSet = ForeignKey(
        FeatureSet, PROTECT, related_name="links_artifact"
    )
    slot: str | None = CharField(max_length=40, null=True)
    feature_ref_is_semantic: bool | None = BooleanField(
        null=True
    )  # like Feature name or Gene symbol or CellMarker name

    class Meta:
        unique_together = ("artifact", "featureset")


class CollectionArtifact(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_artifact"
    )
    artifact: Artifact = ForeignKey(Artifact, PROTECT, related_name="links_collection")

    class Meta:
        unique_together = ("collection", "artifact")


class ArtifactULabel(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_artifactulabel", default=None
    )
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        # can have the same label linked to the same artifact if the feature is
        # different
        unique_together = ("artifact", "ulabel", "feature")


class TransformULabel(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_ulabel")
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_transform")

    class Meta:
        unique_together = ("transform", "ulabel")


class CollectionULabel(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_ulabel"
    )
    ulabel: ULabel = ForeignKey(ULabel, PROTECT, related_name="links_collection")
    feature: Feature | None = ForeignKey(
        Feature, PROTECT, null=True, related_name="links_collectionulabel", default=None
    )
    label_ref_is_name: bool | None = BooleanField(null=True)
    feature_ref_is_name: bool | None = BooleanField(null=True)

    class Meta:
        unique_together = ("collection", "ulabel")


class ArtifactFeatureValue(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="+")
    # we follow the lower() case convention rather than snake case for link models
    featurevalue = ForeignKey(FeatureValue, PROTECT, related_name="+")

    class Meta:
        unique_together = ("artifact", "featurevalue")


class RunParamValue(BasicRecord, LinkORM):
    id: int = models.BigAutoField(primary_key=True)
    run: Run = ForeignKey(Run, CASCADE, related_name="+")
    # we follow the lower() case convention rather than snake case for link models
    paramvalue: ParamValue = ForeignKey(ParamValue, PROTECT, related_name="+")
    created_at: datetime = DateTimeField(auto_now_add=True, db_index=True)
    """Time of creation of record."""
    created_by: User = ForeignKey(
        "lamindb.User", PROTECT, default=current_user_id, related_name="+"
    )
    """Creator of record."""

    class Meta:
        unique_together = ("run", "paramvalue")


class ArtifactParamValue(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="+")
    # we follow the lower() case convention rather than snake case for link models
    paramvalue: ParamValue = ForeignKey(ParamValue, PROTECT, related_name="+")

    class Meta:
        unique_together = ("artifact", "paramvalue")


# -------------------------------------------------------------------------------------
# Link models for project management


class ArtifactProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_artifactproject",
    )
    label_ref_is_name: bool | None = BooleanField(null=True, default=None)
    feature_ref_is_name: bool | None = BooleanField(null=True, default=None)

    class Meta:
        # can have the same label linked to the same artifact if the feature is different
        unique_together = ("artifact", "project", "feature")


class TransformProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(Transform, CASCADE, related_name="links_project")
    project: Project = ForeignKey(Project, PROTECT, related_name="links_transform")

    class Meta:
        unique_together = ("transform", "project")


class CollectionProject(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_project"
    )
    project: Project = ForeignKey(Project, PROTECT, related_name="links_collection")

    class Meta:
        unique_together = ("collection", "project")


class ArtifactReference(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    artifact: Artifact = ForeignKey(Artifact, CASCADE, related_name="links_reference")
    reference: Reference = ForeignKey(Reference, PROTECT, related_name="links_artifact")
    feature: Feature | None = ForeignKey(
        Feature,
        PROTECT,
        null=True,
        default=None,
        related_name="links_artifactreference",
    )
    label_ref_is_name: bool | None = BooleanField(null=True, default=None)
    feature_ref_is_name: bool | None = BooleanField(null=True, default=None)

    class Meta:
        # can have the same label linked to the same artifact if the feature is different
        unique_together = ("artifact", "reference", "feature")


class TransformReference(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    transform: Transform = ForeignKey(
        Transform, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_transform"
    )

    class Meta:
        unique_together = ("transform", "reference")


class CollectionReference(BasicRecord, LinkORM, TracksRun):
    id: int = models.BigAutoField(primary_key=True)
    collection: Collection = ForeignKey(
        Collection, CASCADE, related_name="links_reference"
    )
    reference: Reference = ForeignKey(
        Reference, PROTECT, related_name="links_collection"
    )

    class Meta:
        unique_together = ("collection", "reference")


# class Migration(Record):
#     app = CharField(max_length=255)
#     name = CharField(max_length=255)
#     applied: datetime = DateTimeField()

#     class Meta:
#         db_table = "django_migrations"
#         managed = False


# -------------------------------------------------------------------------------------
# Low-level logic needed in lamindb-setup

# Below is needed within lnschema-core because lamindb-setup already performs
# some logging


def format_field_value(value: datetime | str | Any) -> Any:
    from datetime import datetime

    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M:%S %Z")

    if isinstance(value, str):
        try:
            value = datetime.fromisoformat(value)
            value = value.strftime("%Y-%m-%d %H:%M:%S %Z")
        except ValueError:
            pass
        return f"'{value}'"
    else:
        return value


class RegistryInfo:
    def __init__(self, registry: Registry):
        self.registry = registry

    def _get_type_for_field(self, field_name: str) -> str:
        field = self.registry._meta.get_field(field_name)
        related_model_name = (
            field.related_model.__name__
            if hasattr(field, "related_model") and field.related_model
            else None
        )
        return related_model_name if related_model_name else field.get_internal_type()

    def _get_base_class_fields(self) -> list[str]:
        return [
            field.name
            for base in self.registry.__bases__
            if hasattr(base, "_meta")
            for field in base._meta.get_fields()
        ]

    def _reorder_fields_by_class(self, fields_to_order: list[Field]) -> list[Field]:
        """Reorders the fields so that base class fields come last."""
        non_base_class_fields = [
            field
            for field in fields_to_order
            if field.name not in self._get_base_class_fields()
        ]
        found_base_class_fields = [
            field
            for field in fields_to_order
            if field.name in self._get_base_class_fields()
        ]
        return non_base_class_fields + found_base_class_fields

    def get_simple_fields(self, return_str: bool = False) -> Any:
        simple_fields = [
            field
            for field in self.registry._meta.get_fields()
            if not (
                isinstance(field, ManyToOneRel)
                or isinstance(field, ManyToManyRel)
                or isinstance(field, ManyToManyField)
                or isinstance(field, ForeignKey)
                or field.name.startswith("_")
                or field.name == "id"
            )
        ]
        simple_fields = self._reorder_fields_by_class(simple_fields)
        if not return_str:
            return simple_fields
        else:
            repr_str = f"  {colors.italic('Simple fields')}\n"
            if simple_fields:
                repr_str += "".join(
                    [
                        f"    .{field_name.name}: {self._get_type_for_field(field_name.name)}\n"
                        for field_name in simple_fields
                    ]
                )
            return repr_str

    def get_relational_fields(self, return_str: bool = False):
        # we ignore ManyToOneRel because it leads to so much clutter in the API
        # also note that our general guideline is to have related_name="+"
        # for ForeignKey fields
        relational_fields = (ManyToOneRel, ManyToManyRel, ManyToManyField, ForeignKey)

        class_specific_relational_fields = [
            field
            for field in self.registry._meta.fields + self.registry._meta.many_to_many
            if isinstance(field, relational_fields)
            and not field.name.startswith(("links_", "_"))
        ]

        non_class_specific_relational_fields = [
            field
            for field in self.registry._meta.get_fields()
            if isinstance(field, relational_fields)
            and not field.name.startswith(("links_", "_"))
        ]
        non_class_specific_relational_fields = self._reorder_fields_by_class(
            non_class_specific_relational_fields
        )

        # Ensure that class specific fields (e.g. Artifact) come before non-class specific fields (e.g. collection)
        filtered_non_class_specific = [
            field
            for field in non_class_specific_relational_fields
            if field not in class_specific_relational_fields
        ]
        ordered_relational_fields = (
            class_specific_relational_fields + filtered_non_class_specific
        )

        core_module_fields = []
        external_modules_fields = []
        for field in ordered_relational_fields:
            field_name = repr(field).split(": ")[1][:-1]
            if field_name.count(".") == 1 and "lamindb" not in field_name:
                external_modules_fields.append(field)
            else:
                core_module_fields.append(field)

        def _get_related_field_type(field) -> str:
            field_type = (
                field.related_model.__get_name_with_module__()
                .replace(
                    "Artifact", ""
                )  # some fields have an unnecessary 'Artifact' in their name
                .replace(
                    "Collection", ""
                )  # some fields have an unnecessary 'Collection' in their name
            )
            return (
                self._get_type_for_field(field.name)
                if not field_type.strip()
                else field_type
            )

        core_module_fields_formatted = [
            f"    .{field.name}: {_get_related_field_type(field)}\n"
            for field in core_module_fields
        ]
        external_modules_fields_formatted = [
            f"    .{field.name}: {_get_related_field_type(field)}\n"
            for field in external_modules_fields
        ]

        if not return_str:
            external_modules_fields_by_modules = defaultdict(list)
            for field_str, field in zip(
                external_modules_fields_formatted, external_modules_fields
            ):
                field_type = field_str.split(":")[1].split()[0]
                module_name = field_type.split(".")[0]
                external_modules_fields_by_modules[module_name].append(field)
            return core_module_fields, external_modules_fields_by_modules
        else:
            repr_str = ""

            # Non-external relational fields
            if core_module_fields:
                repr_str += f"  {colors.italic('Relational fields')}\n"
                repr_str += "".join(core_module_fields_formatted)

            # External relational fields
            external_modules = set()
            for field in external_modules_fields_formatted:
                field_type = field.split(":")[1].split()[0]
                external_modules.add(field_type.split(".")[0])

            if external_modules:
                # We want Bionty to show up before other modules
                external_modules = (
                    ["bionty"] + sorted(external_modules - {"bionty"})  # type: ignore
                    if "bionty" in external_modules
                    else sorted(external_modules)
                )
                for ext_module in external_modules:
                    ext_module_fields = [
                        field
                        for field in external_modules_fields_formatted
                        if ext_module in field
                    ]

                    if ext_module_fields:
                        repr_str += (
                            f"  {colors.italic(f'{ext_module.capitalize()} fields')}\n"
                        )
                        repr_str += "".join(ext_module_fields)

            return repr_str


def registry_repr(cls):
    """Shows fields."""
    repr_str = f"{colors.green(cls.__name__)}\n"
    info = RegistryInfo(cls)
    repr_str += info.get_simple_fields(return_str=True)
    repr_str += info.get_relational_fields(return_str=True)
    repr_str = repr_str.rstrip("\n")
    return repr_str


def record_repr(
    self: Record, include_foreign_keys: bool = True, exclude_field_names=None
) -> str:
    if exclude_field_names is None:
        exclude_field_names = ["id", "updated_at", "source_code"]
    field_names = [
        field.name
        for field in self._meta.fields
        if (not isinstance(field, ForeignKey) and field.name not in exclude_field_names)
    ]
    if include_foreign_keys:
        field_names += [
            f"{field.name}_id"
            for field in self._meta.fields
            if isinstance(field, ForeignKey)
        ]
    if "created_at" in field_names:
        field_names.remove("created_at")
        field_names.append("created_at")
    if field_names[0] != "uid" and "uid" in field_names:
        field_names.remove("uid")
        field_names.insert(0, "uid")
    fields_str = {}
    for k in field_names:
        if not k.startswith("_") and hasattr(self, k):
            value = getattr(self, k)
            # Force strip the time component of the version
            if k == "version" and value:
                fields_str[k] = f"'{str(value).split()[0]}'"
            else:
                fields_str[k] = format_field_value(value)
    fields_joined_str = ", ".join(
        [f"{k}={fields_str[k]}" for k in fields_str if fields_str[k] is not None]
    )
    return f"{self.__class__.__name__}({fields_joined_str})"


# below is code to further format the repr of a record
#
# def format_repr(
#     record: Record, exclude_field_names: str | list[str] | None = None
# ) -> str:
#     if isinstance(exclude_field_names, str):
#         exclude_field_names = [exclude_field_names]
#     exclude_field_names_init = ["id", "created_at", "updated_at"]
#     if exclude_field_names is not None:
#         exclude_field_names_init += exclude_field_names
#     return record.__repr__(
#         include_foreign_keys=False, exclude_field_names=exclude_field_names_init
#     )


Record.__repr__ = record_repr  # type: ignore
Record.__str__ = record_repr  # type: ignore


def deferred_attribute__repr__(self):
    return f"FieldAttr({self.field.model.__name__}.{self.field.name})"


FieldAttr.__repr__ = deferred_attribute__repr__  # type: ignore
# backward compatibility
CanValidate = CanCurate
