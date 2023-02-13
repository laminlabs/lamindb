# it'd be nice to use properties, but class properties are
# a difficult thing
# all hacks via a global instance or a metaclass lead to new hacks
# in the docs, so, we'll go with the below
class settings:
    """Settings.

    This is a static class to manage post-setup settings.

    See `lamindb.settings <https://lamin.ai/docs/db/lamindb.settings>`__.

    For setup-related settings, see
    `lndb.settings <https://lamin.ai/docs/lndb/lndb.settings>`__.
    """

    error_on_dobject_hash_exists: bool = True
    """Upon ingestion, error if a dobject hash equals an existing hash in the DB.

    FAQ: :doc:`/faq/ingest-same-file-twice`.
    """
    track_run_inputs_upon_load: bool = False
    """Upon load, add loaded dobjects as the input of the current notebook run.

    FAQ: :doc:`/faq/track-runin`.
    """
