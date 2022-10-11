from typing import Dict

from ...dev.file import load_to_memory
from ...schema import wetlab
from .._link import link
from .._query import query


class LinkStaged:
    """Link db entries to the dobject, accessible via `Staged().link`.

    Args:
        ingest: an :class:`~lamindb.db.ingest.Staged` instance.
    """

    def __init__(self, ingest) -> None:
        self._ingest = ingest
        # TODO: make sure dtransform is unique
        self._entries: Dict = {}  # staged entries to be inserted

    @property
    def linked_entries(self) -> dict:
        """Linked db entries of the dobject."""
        return self._entries

    def add_entry(self, table_name: str, entry) -> None:
        """Add an entry to ._entries."""
        if table_name not in self._entries:
            self._entries[table_name] = []
        self._entries[table_name].append(entry)

    def features(self, feature_model, *, featureset_name: str = None) -> None:
        """Link dobject to features.

        Args:
            feature_model: a feature model instance
            featureset_name: name of the featureset

        Returns:
            writes to `Staged.feature_model`
        """
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

    def biometa(self, biometa: wetlab.biometa) -> None:
        """Link dobject to a biometa.

        Args:
            biometa: a wetlab.biometa entry
        """
        result = query.biometa(id=biometa.id).one_or_none()  # type: ignore
        if result is None:
            self.add_entry("biometa", biometa)
        dobject_id = self._ingest.dobject.id

        # create an entry in the link table
        link_entry = query.dobject_biometa(  # type: ignore
            dobject_id=dobject_id, biometa_id=biometa.id
        ).one_or_none()
        if link_entry is None:
            link_entry = wetlab.dobject_biometa(
                dobject_id=dobject_id, biometa_id=biometa.id
            )
            self.add_entry("dobject_biometa", link_entry)
