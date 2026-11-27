from laminprofiler import laminprofiler


@laminprofiler("RMlmBHcPqL3H0M0K", fail_seconds=1.5)
def profile_import_lamindb_core_storage():
    import lamindb.core.storage  # noqa: F401
