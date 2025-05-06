"""Examples.

.. autosummary::
   :toctree: .

   ingest_mini_immuno_datasets
   schemas

"""

from . import schemas


def ingest_mini_immuno_datasets():
    """Ingest mini immuno datasets.

    .. literalinclude:: scripts/ingest_mini_immuno_datasets.py
        :language: python
    """
    import sys
    from pathlib import Path

    docs_path = Path(__file__).parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    import ingest_mini_immuno_datasets  # noqa
