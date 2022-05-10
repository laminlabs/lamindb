class Configure:
    """Configuration."""

    def __init__(self):
        try:
            from lamindb._configuration import storage_root
        except ImportError:
            raise RuntimeError("Please run: lamin configure")
        self._storage_root = storage_root

    @property
    def storage_root(self):
        return self._storage_root
