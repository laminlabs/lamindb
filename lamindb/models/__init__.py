"""Models library.

.. autosummary::
   :toctree: .

   BasicRecord
   Record
   Registry
   Space
   QuerySet
   QueryManager
   RecordList
   FeatureManager
   ParamManager
   LabelManager
   IsVersioned
   CanCurate
   HasParents
   TracksRun
   TracksUpdates
   ParamValue
   FeatureValue
   InspectResult
   ValidateFields

"""

# ruff: noqa: I001
from lamin_utils._inspect import InspectResult
from ._is_versioned import IsVersioned
from .can_curate import CanCurate
from .record import (
    BasicRecord,
    Record,
    Registry,
    Space,
    ValidateFields,
    format_field_value,
    record_repr,
    LinkORM,
)
from .core import Storage
from .transform import Transform
from .run import Run, TracksRun, TracksUpdates, current_run, Param, ParamValue, User
from .feature import Feature, FeatureValue
from .schema import Schema
from .ulabel import ULabel

# should come last as it needs everything else
from .artifact import Artifact
from ._feature_manager import FeatureManager
from .run import ParamManager
from ._label_manager import LabelManager
from .collection import Collection, CollectionArtifact
from .project import Person, Project, Reference
from .flextable import FlexTable, RunData
from .query_manager import QueryManager
from .query_set import QuerySet, RecordList
from .has_parents import HasParents
from datetime import datetime as _datetime

FeatureSet = Schema  # backward compat
