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
.. autoclass:: ArtifactJsonValue
.. autoclass:: ArtifactProject
.. autoclass:: ArtifactReference
.. autoclass:: ArtifactRun
.. autoclass:: ArtifactSchema
.. autoclass:: ArtifactULabel
.. autoclass:: ArtifactUser
.. autoclass:: BlockProject
.. autoclass:: CollectionArtifact
.. autoclass:: CollectionProject
.. autoclass:: CollectionReference
.. autoclass:: CollectionULabel
.. autoclass:: FeatureProject
.. autoclass:: RunJsonValue
.. autoclass:: RunProject
.. autoclass:: RunULabel
.. autoclass:: SchemaComponent
.. autoclass:: SchemaFeature
.. autoclass:: SchemaProject
.. autoclass:: TransformProject
.. autoclass:: TransformReference
.. autoclass:: TransformTransform
.. autoclass:: TransformULabel
.. autoclass:: ULabelProject

Record values
-------------

Record values work almost exactly like artifact and run annotations,
with the exception that JSON values are stored in `RecordJson` on a per-record basis
and not in `JsonValue`.

.. autoclass:: ArtifactRecord
.. autoclass:: CollectionRecord
.. autoclass:: ProjectRecord
.. autoclass:: RecordArtifact
.. autoclass:: RecordCollection
.. autoclass:: RecordJson
.. autoclass:: RecordProject
.. autoclass:: RecordRecord
.. autoclass:: RecordReference
.. autoclass:: RecordRun
.. autoclass:: RecordTransform
.. autoclass:: RecordULabel
.. autoclass:: RecordUser
.. autoclass:: ReferenceRecord
.. autoclass:: RunRecord
.. autoclass:: TransformRecord

Blocks
------

.. autoclass:: Block
.. autoclass:: ArtifactBlock
.. autoclass:: BranchBlock
.. autoclass:: CollectionBlock
.. autoclass:: FeatureBlock
.. autoclass:: ProjectBlock
.. autoclass:: RecordBlock
.. autoclass:: RunBlock
.. autoclass:: SchemaBlock
.. autoclass:: SpaceBlock
.. autoclass:: TransformBlock
.. autoclass:: ULabelBlock

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
from .transform import Transform, TransformTransform
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
    ArtifactReference,
    BlockProject,
    CollectionProject,
    CollectionReference,
    FeatureProject,
    ProjectRecord,
    RecordProject,
    RecordReference,
    ReferenceRecord,
    RunProject,
    SchemaProject,
    TransformProject,
    TransformReference,
    ULabelProject,
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
    ArtifactRecord,
    CollectionRecord,
    RecordArtifact,
    RecordCollection,
    RecordJson,
    RecordRecord,
    RecordRun,
    RecordTransform,
    RecordULabel,
    RecordUser,
    RunRecord,
    TransformRecord,
)
from .block import (
    Block,
    ArtifactBlock,
    BranchBlock,
    CollectionBlock,
    FeatureBlock,
    ProjectBlock,
    RecordBlock,
    RunBlock,
    SchemaBlock,
    SpaceBlock,
    TransformBlock,
    ULabelBlock,
)

FeatureValue = JsonValue  # backward compatibility
