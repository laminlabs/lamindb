from pathlib import Path


# We could use pydantic BaseSettings
# I just don't understand how we would then raise an informative error
# guiding the user to run the configuration dialogue
class settings:
    @classmethod
    @property
    def storage_root(cls) -> Path:
        """Root directory of storage (local or cloud)."""
        try:
            from lamindb._configuration import storage_root
        except ImportError:
            raise RuntimeError(
                "Please configure lamindb in the CLI. Run: lamindb configure"
            )
        return Path(storage_root)
