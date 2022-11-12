"""Manage Jupyter notebooks.

You'll typically only need one function in this module: `header()` initializes
and prints metadata for a Jupyter notebook.

For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.

.. autosummary::
   :toctree: .

   header
   publish
   run
"""
from ._nb import _run, header, publish  # noqa

run = _run
"""Current notebook run."""
