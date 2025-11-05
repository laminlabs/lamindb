"""Models library.

Feature and label managers
--------------------------

.. autoclass:: FeatureManager
.. autoclass:: LabelManager

Registry base classes
---------------------

.. autoclass:: BaseSQLRecord
.. autoclass:: SQLRecord
.. autoclass:: Registry

Mixins for registries
---------------------

.. autoclass:: IsVersioned
.. autoclass:: CanCurate
.. autoclass:: HasParents
.. autoclass:: TracksRun
.. autoclass:: TracksUpdates

Query sets & managers
---------------------

.. autoclass:: BasicQuerySet
.. autoclass:: QuerySet
.. autoclass:: ArtifactSet
.. autoclass:: QueryManager

Storage of feature values
-------------------------

.. autoclass:: FeatureValue

Utility classes
---------------

.. autoclass:: LazyArtifact
.. autoclass:: SQLRecordList
.. autoclass:: InspectResult
.. autoclass:: ValidateFields
.. autoclass:: SchemaOptionals

"""

# ruff: noqa: I001
from lamin_utils._inspect import InspectResult
from ._is_versioned import IsVersioned
from .can_curate import CanCurate
from .sqlrecord import (
    BaseSQLRecord,
    SQLRecord,
    Registry,
    Space,
    Branch,
    Migration,
    ValidateFields,
    format_field_value,
    record_repr,
    IsLink,
)
from .storage import Storage
from .transform import Transform
from .run import Run, TracksRun, TracksUpdates, current_run, User
from .feature import Feature, FeatureValue
from .schema import Schema
from .ulabel import ULabel

# should come last as it needs everything else
from .artifact import Artifact, LazyArtifact
from ._feature_manager import FeatureManager
from ._label_manager import LabelManager
from .collection import Collection, CollectionArtifact
from .project import Project, Reference
from .query_manager import QueryManager
from .query_set import BasicQuerySet, QuerySet, SQLRecordList
from .artifact_set import ArtifactSet
from .has_parents import HasParents
from datetime import datetime as _datetime

FeatureSet = Schema  # backward compat

# link models
from .artifact import ArtifactFeatureValue
from .project import (
    ArtifactProject,
    TransformProject,
    CollectionProject,
    ULabelProject,
    FeatureProject,
    SchemaProject,
    ArtifactReference,
    CollectionReference,
    RunProject,
    RecordProject,
    RecordReference,
    ReferenceRecord,
    ProjectRecord,
)
from .run import RunFeatureValue
from .schema import (
    SchemaFeature,
    ArtifactSchema,
    SchemaComponent,
    SchemaOptionals,
)
from .ulabel import ArtifactULabel, TransformULabel, RunULabel, CollectionULabel

from .record import (
    Record,
    RecordJson,
    RecordRecord,
    RecordULabel,
    RecordRun,
    RunRecord,
    RecordUser,
    RecordArtifact,
    ArtifactRecord,
)
from .block import (
    RootBlock,
    ArtifactBlock,
    TransformBlock,
    RecordBlock,
    CollectionBlock,
    RunBlock,
    SchemaBlock,
    ProjectBlock,
    BranchBlock,
    SpaceBlock,
)

LinkORM = IsLink  # backward compat
ParamValue = FeatureValue  # backward compat
ArtifactParamValue = ArtifactFeatureValue  # backward compat
RunParamValue = RunFeatureValue  # backward compat
Param = Feature  # backward compat
BasicRecord = BaseSQLRecord  # backward compat
