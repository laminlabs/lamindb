from laminprofiler import laminprofiler


@laminprofiler("avOG1AOJCSq7nooG", fail_seconds=3.5)
def profile_import_lamindb_and_connect():
    import lamindb as ln

    # Log the profiling run on laminlabs/lamindb-benchmarks, so connect elsewhere.
    ln.connect("laminlabs/lamin-site-assets")
