"""Manage Jupyter notebooks.

.. autosummary::
   :toctree: .

   run
   header
   publish

For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.
"""
from ._nb import _run as run  # noqa
from ._nb import header, publish  # noqa
