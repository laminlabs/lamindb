"""Manage Jupyter notebooks.

Two functions:

.. autosummary::
   :toctree: .

   header
   publish

Two records:

.. autosummary::
   :toctree: .

   jupynb
   run

For more background, see `nbproject <https://lamin.ai/docs/nbproject>`__.

"""
from ._nb import _jupynb, _run, header, publish  # noqa

jupynb = _jupynb
"""Current notebook record."""
run = _run
"""Current run record."""
