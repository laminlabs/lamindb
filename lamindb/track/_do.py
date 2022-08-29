import lamindb as lndb


def do():
    """Load data access log."""
    lndb.db.load("usage")
