from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, final

import lamindb_setup as ln_setup
from lamin_utils import logger
from lamindb_setup import settings as setup_settings

if TYPE_CHECKING:
    import pandas as pd
    from bionty.models import (
        CellLine,
        CellMarker,
        CellType,
        DevelopmentalStage,
        Disease,
        Ethnicity,
        ExperimentalFactor,
        Gene,
        Organism,
        Pathway,
        Phenotype,
        Protein,
        Tissue,
    )
    from pertdb.models import (
        Biologic,
        CombinationPerturbation,
        Compound,
        CompoundPerturbation,
        EnvironmentalPerturbation,
        GeneticPerturbation,
        PerturbationTarget,
    )

    from lamindb.models import (
        Artifact,
        Branch,
        Collection,
        Feature,
        Project,
        QuerySet,
        Record,
        Reference,
        Run,
        Schema,
        Space,
        Storage,
        Transform,
        ULabel,
        User,
    )

    from .sqlrecord import SQLRecord


@final
class NonInstantiableQuerySet:
    """Wrapper around QuerySet that prevents instantiation while preserving query methods."""

    def __init__(self, qs: QuerySet, registry_name: str):
        self._qs = qs
        self._name = registry_name

    def __repr__(self) -> str:
        return f"<QuerySet [{self._name}]>"

    def __call__(self, *args, **kwargs):
        raise TypeError(
            f"Cannot instantiate {self._name} from DB. "
            f"Use {self._name}.filter(), {self._name}.get(), etc. to query records."
        )

    def __getattr__(self, attr):
        return getattr(self._qs, attr)


class ModuleNamespace:
    """Namespace for accessing registries from a specific schema module.

    Args:
        query_db: Parent DB instance.
        module_name: Name of the schema module (e.g., 'bionty', 'pertdb').
    """

    def __init__(self, query_db: DB, module_name: str):
        self._query_db = query_db
        self._module_name = module_name
        self._cache: dict[str, NonInstantiableQuerySet] = {}

    def __getattr__(self, name: str) -> NonInstantiableQuerySet:
        """Access a registry class from this schema module.

        Args:
            name: Registry class name (e.g., 'Gene', 'CellType').

        Returns:
            QuerySet for the specified registry scoped to the parent instance.
        """
        if name in self._cache:
            return self._cache[name]

        try:
            schema_module = import_module(self._module_name)
            if hasattr(schema_module, name):
                model_class = getattr(schema_module, name)
                queryset = model_class.connect(self._query_db._instance)
                wrapped = NonInstantiableQuerySet(queryset, name)
                self._cache[name] = wrapped
                return wrapped
        except (ImportError, AttributeError):
            pass

        raise AttributeError(
            f"Registry '{name}' not found in lamindb. Use .bt.{name} or .pertdb.{name} for schema-specific registries."
        )

    def __dir__(self) -> list[str]:
        """Return list of available registries in this schema module."""
        base_attrs = [attr for attr in object.__dir__(self) if not attr.startswith("_")]
        try:
            schema_module = import_module(self._module_name)
            if hasattr(schema_module, "__all__"):
                registries = set()
                for class_name in schema_module.__all__:
                    model_class = getattr(schema_module, class_name, None)
                    if model_class and hasattr(model_class, "connect"):
                        registries.add(class_name)
                return sorted(set(base_attrs) | registries)
        except ImportError:
            pass
        return base_attrs


class BiontyDB(ModuleNamespace):
    """Namespace for Bionty registries (Gene, CellType, Disease, etc.)."""

    Gene: QuerySet[Gene]  # type: ignore[type-arg]
    Protein: QuerySet[Protein]  # type: ignore[type-arg]
    CellType: QuerySet[CellType]  # type: ignore[type-arg]
    Disease: QuerySet[Disease]  # type: ignore[type-arg]
    Phenotype: QuerySet[Phenotype]  # type: ignore[type-arg]
    Pathway: QuerySet[Pathway]  # type: ignore[type-arg]
    Tissue: QuerySet[Tissue]  # type: ignore[type-arg]
    CellLine: QuerySet[CellLine]  # type: ignore[type-arg]
    CellMarker: QuerySet[CellMarker]  # type: ignore[type-arg]
    Organism: QuerySet[Organism]  # type: ignore[type-arg]
    ExperimentalFactor: QuerySet[ExperimentalFactor]  # type: ignore[type-arg]
    DevelopmentalStage: QuerySet[DevelopmentalStage]  # type: ignore[type-arg]
    Ethnicity: QuerySet[Ethnicity]  # type: ignore[type-arg]


class PertdbDB(ModuleNamespace):
    """Namespace for `PertDB` registries (Biologic, Compound, etc.)."""

    Biologic: QuerySet[Biologic]  # type: ignore[type-arg]
    Compound: QuerySet[Compound]  # type: ignore[type-arg]
    CompoundPerturbation: QuerySet[CompoundPerturbation]  # type: ignore[type-arg]
    GeneticPerturbation: QuerySet[GeneticPerturbation]  # type: ignore[type-arg]
    EnvironmentalPerturbation: QuerySet[EnvironmentalPerturbation]  # type: ignore[type-arg]
    CombinationPerturbation: QuerySet[CombinationPerturbation]  # type: ignore[type-arg]
    PerturbationTarget: QuerySet[PerturbationTarget]  # type: ignore[type-arg]


class DB:
    """Query any registry of any instance.

    Args:
        instance: Instance identifier in format "account/instance".

    Examples
    --------

    Query objects from an instance::

        db = ln.DB("laminlabs/cellxgene")

    Query artifacts and filter by `suffix`::

        db.Artifact.filter(suffix=".h5ad").to_dataframe()

    Get a single artifact by uid::

        artifact = db.Artifact.get("abcDEF123456")

    Query records and filter by name::

        db.Record.filter(name__startswith="sample").to_dataframe()

    Get a cell type object::

        t_cell = db.bionty.CellType.get(name="T cell")

    Create a lookup object to auto-complete all cell types in the database::

        cell_types = db.bionty.CellType.lookup()

    Return a `DataFrame` with additional info::

        db.Artifact.filter(
            suffix=".h5ad",
            description__contains="immune",
            size__gt=1e9,  # size > 1GB
            cell_types__name__in=["B cell", "T cell"],
        ).order_by("created_at").to_dataframe(
            include=["cell_types__name", "created_by__handle"]  # include additional info
        ).head()
    """

    Artifact: QuerySet[Artifact]  # type: ignore[type-arg]
    Collection: QuerySet[Collection]  # type: ignore[type-arg]
    Transform: QuerySet[Transform]  # type: ignore[type-arg]
    Run: QuerySet[Run]  # type: ignore[type-arg]
    User: QuerySet[User]  # type: ignore[type-arg]
    Storage: QuerySet[Storage]  # type: ignore[type-arg]
    Feature: QuerySet[Feature]  # type: ignore[type-arg]
    ULabel: QuerySet[ULabel]  # type: ignore[type-arg]
    Record: QuerySet[Record]  # type: ignore[type-arg]
    Schema: QuerySet[Schema]  # type: ignore[type-arg]
    Project: QuerySet[Project]  # type: ignore[type-arg]
    Reference: QuerySet[Reference]  # type: ignore[type-arg]
    Branch: QuerySet[Branch]  # type: ignore[type-arg]
    Space: QuerySet[Space]  # type: ignore[type-arg]

    bionty: BiontyDB
    pertdb: PertdbDB

    def __init__(self, instance: str):
        self._instance = instance
        self._cache: dict[str, NonInstantiableQuerySet | BiontyDB | PertdbDB] = {}
        self._available_registries: set[str] | None = None

        owner, instance_name = (
            ln_setup._connect_instance.get_owner_name_from_identifier(instance)
        )
        instance_info = ln_setup._connect_instance._connect_instance(
            owner=owner,
            name=instance_name,
            allow_sqlite_clone_fallback=True,
        )
        self._instance_info = instance_info
        self._modules = ["lamindb"] + list(instance_info.modules)
        warning = ln_setup.core.django._warn_module_mismatch(
            target_apps={"lamindb"} | instance_info.modules,
            # Read-only DB querying should only warn when instance modules are missing
            # from the local environment, not when local modules are additional.
            current_apps={"lamindb"} | (setup_settings.modules & instance_info.modules),
        )
        if warning is not None:
            logger.warning(warning)

    def __getattr__(self, name: str) -> NonInstantiableQuerySet | BiontyDB | PertdbDB:
        """Access a registry class or schema namespace for this database instance.

        Args:
            name: Registry class name (e.g., 'Artifact', 'Collection') or schema namespace ('bionty', 'pertdb').

        Returns:
            QuerySet for the specified registry or schema namespace scoped to this instance.
        """
        if name in self._cache:
            return self._cache[name]

        if name == "bionty":
            if "bionty" not in self._modules:
                raise AttributeError(
                    f"Schema 'bionty' not available in instance '{self._instance}'."
                )
            if "bionty" not in self._cache:
                namespace = BiontyDB(self, "bionty")
                self._cache["bionty"] = namespace
            return self._cache["bionty"]

        if name == "pertdb":
            if "pertdb" not in self._modules:
                raise AttributeError(
                    f"Schema 'pertdb' not available in instance '{self._instance}'."
                )
            if "pertdb" not in self._cache:
                namespace = PertdbDB(self, "pertdb")  # type: ignore
                self._cache["pertdb"] = namespace
            return self._cache["pertdb"]

        lamindb_module = import_module("lamindb")
        if hasattr(lamindb_module, name):
            model_class = getattr(lamindb_module, name)
            queryset = model_class.connect(
                self._instance, _instance_info=self._instance_info
            )
            wrapped = NonInstantiableQuerySet(queryset, name)
            self._cache[name] = wrapped
            return wrapped
        else:
            raise AttributeError(
                f"Registry '{name}' not found in lamindb core registries. Use .bionty.{name} or .pertdb.{name} for schema-specific registries."
            )

    def view(
        self,
        *,
        limit: int = 7,
        modules: str | None = None,
        registries: list[str] | None = None,
        df: pd.DataFrame | None = None,
    ) -> None:
        """View metadata for this database instance."""
        from ._view import _view

        def get_queryable(registry: type[SQLRecord], module_name: str):
            if module_name == "core":
                return getattr(self, registry.__name__)
            return getattr(getattr(self, module_name), registry.__name__)

        return _view(
            limit=limit,
            modules=modules,
            registries=registries,
            df=df,
            default_modules=["core"] + list(self._instance_info.modules),
            get_schema_module=lambda module_name: import_module(
                "lamindb" if module_name == "core" else module_name
            ),
            get_queryable=get_queryable,
        )

    def __repr__(self) -> str:
        return f"DB('{self._instance}')"

    def __dir__(self) -> list[str]:
        """Return list of available registries and schema namespaces."""
        base_attrs = [attr for attr in super().__dir__() if not attr.startswith("_")]

        lamindb_registries = set()
        try:
            lamindb_module = import_module("lamindb")
            if hasattr(lamindb_module, "__all__"):
                for class_name in lamindb_module.__all__:
                    model_class = getattr(lamindb_module, class_name, None)
                    if model_class and hasattr(model_class, "connect"):
                        lamindb_registries.add(class_name)
        except ImportError:
            pass

        module_namespaces = set()
        if "bionty" in self._modules:
            module_namespaces.add("bionty")
        if "pertdb" in self._modules:
            module_namespaces.add("pertdb")

        return sorted(set(base_attrs) | lamindb_registries | module_namespaces)
