from laminprofiler import laminprofiler


@laminprofiler("iJEIHhE5OVrz7qmP", total_seconds=1.5)
def profile_import_lamindb():
    import lamindb  # noqa: F401
