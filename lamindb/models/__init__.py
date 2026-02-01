"""Auxiliary models & database library.

Registry basics
---------------

.. autoclass:: BaseSQLRecord
.. autoclass:: SQLRecord
.. autoclass:: Registry
.. autoclass:: QuerySet

Mixins for registries
---------------------

.. autoclass:: IsVersioned
.. autoclass:: HasType
.. autoclass:: HasParents
.. autoclass:: CanCurate
.. autoclass:: TracksRun
.. autoclass:: TracksUpdates

Managers
--------

.. autoclass:: FeatureManager
.. autoclass:: LabelManager
.. autoclass:: QueryManager
.. autoclass:: RelatedManager

Artifact & run annotations
--------------------------

Artifact & run annotations can be conditioned on features.
Besides linking categorical data, you can also link simple data types
by virtue of the `JsonValue` model.

.. autoclass:: JsonValue

Below follow the underlying link models for annotations.

.. autoclass:: ArtifactArtifact
.. autoclass:: ArtifactReference
.. autoclass:: ArtifactProject
.. autoclass:: ArtifactTransform
.. autoclass:: ArtifactRun
.. autoclass:: ArtifactCollection
.. autoclass:: ArtifactSchema
.. autoclass:: ArtifactUser
.. autoclass:: ArtifactULabel
.. autoclass:: ArtifactRecord
.. autoclass:: ArtifactReference
.. autoclass:: ArtifactProject
.. autoclass:: ArtifactTransform
.. autoclass:: ArtifactRun
.. autoclass:: ArtifactCollection
.. autoclass:: ArtifactSchema
.. autoclass:: ArtifactUser
.. autoclass:: ArtifactULabel
.. autoclass:: RunProject
.. autoclass:: RunReference
.. autoclass:: RunTransform
.. autoclass:: RunRun
.. autoclass:: RunCollection
.. autoclass:: RunSchema
.. autoclass:: RunUser
.. autoclass:: RunULabel
.. autoclass:: RunRecord
.. autoclass:: RunReference

Record values
-------------

Record values work almost exactly like artifact and run annotations,
with the exception of JSON values are being stored in `RecordJson` on a per-record basis
and not in `JsonValue`.

.. autoclass:: RecordRecord
.. autoclass:: RecordJson
.. autoclass:: RecordULabel
.. autoclass:: RecordRun
.. autoclass:: RecordArtifact
.. autoclass:: RecordReference
.. autoclass:: RecordProject

Blocks
------

.. autoclass:: Block
.. autoclass:: ArtifactBlock
.. autoclass:: TransformBlock
.. autoclass:: RecordBlock
.. autoclass:: CollectionBlock
.. autoclass:: RunBlock
.. autoclass:: SchemaBlock
.. autoclass:: ProjectBlock
.. autoclass:: BranchBlock
.. autoclass:: SpaceBlock
.. autoclass:: RecordUser

Utils
-----

.. autoclass:: LazyArtifact
.. autoclass:: InspectResult
.. autoclass:: ValidateFields
.. autoclass:: SchemaOptionals
.. autoclass:: lamindb.models.query_set.BiontyDB
.. autoclass:: lamindb.models.query_set.PertdbDB

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
    IsLink,
    HasType,
)
from .storage import Storage
from .transform import Transform
from .run import Run, TracksRun, TracksUpdates, current_run, User
from .feature import Feature, JsonValue
from .schema import Schema
from .ulabel import ULabel

# should come last as it needs everything else
from .artifact import Artifact, LazyArtifact
from ._feature_manager import FeatureManager
from ._label_manager import LabelManager
from .collection import Collection, CollectionArtifact
from .project import Project, Reference
from .query_manager import RelatedManager, QueryManager
from .query_set import BasicQuerySet, QuerySet, DB, SQLRecordList
from .artifact_set import ArtifactSet
from .has_parents import HasParents
from datetime import datetime as _datetime

# link models
from .artifact import ArtifactJsonValue, ArtifactArtifact
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
from .run import RunJsonValue
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
    Block,
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

FeatureValue = JsonValue  # backward compatibility
