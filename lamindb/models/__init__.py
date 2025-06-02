"""Models library.

.. autosummary::
   :toctree: .

   BaseSQLRecord
   SQLRecord
   Registry
   BasicQuerySet
   QuerySet
   ArtifactSet
   QueryManager
   SQLRecordList
   FeatureManager
   LabelManager
   IsVersioned
   CanCurate
   HasParents
   TracksRun
   TracksUpdates
   FeatureValue
   InspectResult
   ValidateFields
   SchemaOptionals

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
from .core import Storage
from .transform import Transform
from .run import Run, TracksRun, TracksUpdates, current_run, User
from .feature import Feature, FeatureValue
from .schema import Schema
from .ulabel import ULabel

# should come last as it needs everything else
from .artifact import Artifact
from ._feature_manager import FeatureManager
from ._label_manager import LabelManager
from .collection import Collection, CollectionArtifact
from .project import Person, Project, Reference
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
    SheetProject,
    RunProject,
    RecordProject,
    PersonProject,
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
    Sheet,
    RecordJson,
    RecordRecord,
    RecordULabel,
    RecordRun,
    RecordArtifact,
)


LinkORM = IsLink  # backward compat
ParamValue = FeatureValue  # backward compat
ArtifactParamValue = ArtifactFeatureValue  # backward compat
RunParamValue = RunFeatureValue  # backward compat
Param = Feature  # backward compat
BasicRecord = BaseSQLRecord  # backward compat
