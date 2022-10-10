from typing import Dict

from ...dev.file import load_to_memory
from ...schema import core
from .._link import link
from .._query import query


class LinkIngest:
    """Link db entries to the dobject, accessible via `Ingest().link`.

    Args:
        ingest: an :class:`~lamindb.db.ingest.Ingest` instance.
    """

    def __init__(self, ingest) -> None:
        self._ingest = ingest
        # TODO: need to allow multiple entries from the same table
        self._entries: Dict = {}  # staged entries to be inserted

    @property
    def linked_entries(self) -> dict:
        """Linked db entries of the dobject."""
        return self._entries

    def features(self, feature_model, *, featureset_name: str = None):
        """Link dobject to features."""
        # curate features
        # looks for the id column, if none is found, will assume in the index
        if self._ingest.dmem is None:
            self._ingest._dmem = load_to_memory(self._ingest.data)
        try:
            df = getattr(self._ingest.dmem, "var")  # for AnnData objects
            if callable(df):
                df = self._ingest.dmem
        except AttributeError:
            df = self._ingest.dmem
        # insert feature entries
        # TODO: needs to be staged without inserting here
        self._ingest._feature_model = link.feature_model(
            df=df, feature_model=feature_model, featureset_name=featureset_name
        )

    def pipeline_run(self, pipeline_run: core.pipeline_run) -> core.dtransform:
        """Link dobject to a pipeline run.

        - Create a dtransform entry from pipeline_run (no insert)
        - Link dobject to dtransform
        """
        result = query.pipeline_run(id=pipeline_run.id).one_or_none()  # type: ignore
        if result is None:
            self._entries["pipeline_run"] = pipeline_run
            dtransform = core.dtransform(pipeline_run_id=pipeline_run.id)
        else:
            dtransform = query.dtransform(  # type: ignore
                pipeline_run_id=pipeline_run.id
            ).one()

        self._entries["dtransform"] = dtransform
        self._ingest.dtransform = dtransform

    def jupynb(self, jupynb: core.jupynb) -> core.dtransform:
        """Link dobject to a jupynb.

        - Create a dtransform entry from jupynb (no insert)
        - Link dobject to dtransform
        """
        result = query.jupynb(id=jupynb.id, v=jupynb.v).one_or_none()  # type: ignore
        if result is None:
            self._entries["jupynb"] = jupynb
            dtransform = core.dtransform(jupynb_id=jupynb.id, jupynb_v=jupynb.v)
        else:
            dtransform = query.dtransform(  # type: ignore
                jupynb_id=jupynb.id, jupynb_v=jupynb.v
            ).one()

        self._entries["dtransform"] = dtransform
        self._ingest.dtransform = dtransform
