# it'd be nice to use properties, but class properties are
# a difficult thing
# all hacks via a global instance or a metaclass lead to new hacks
# in the docs, so, we'll go with the below
class settings:
    """Settings.

    This is a static class to manage post-setup settings.

    FAQ: :doc:`/faq/settings`.

    For setup-related settings, see
    `lndb_setup.settings <https://lamin.ai/docs/lndb-setup/lndb_setup.settings>`__.
    """

    error_on_dobject_hash_exists: bool = True
    """Upon ingestion, error if a dobject hash equals an existing hash in the DB.

    FAQ: :doc:`/faq/ingest-same-file-twice`.
    """
