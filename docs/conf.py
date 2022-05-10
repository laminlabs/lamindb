from datetime import datetime
import lamin  # noqa
import sys
from pathlib import Path
from sphinx.application import Sphinx

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]
import lamin  # noqa

for generated in HERE.glob("lamin.*.rst"):
    generated.unlink()

project = "Lamin"
author = "Lamin Dev"
copyright = f"{datetime.now():%Y}, {author}"
release = lamin.__version__.replace(".dirty", "")

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/lamindev/lamin",
    "use_repository_button": True,
}
source_suffix = [".rst", ".ipynb", ".md"]
master_doc = "index"
default_role = "literal"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]
templates_path = ["_templates"]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "IPython.sphinxext.ipython_console_highlighting",  # noqa https://github.com/spatialaudio/nbsphinx/issues/24
    "myst_nb",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]


# Generate the API documentation when building
jupyter_execute_notebooks = "off"
autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
nitpicky = True  # Report broken links
nitpick_ignore = [
    ("py:meth", "pandas.DataFrame.iloc"),
    ("py:meth", "pandas.DataFrame.loc"),
]
suppress_warnings = ["ref.citation"]


def setup(app: Sphinx):
    # Don’t allow broken links. DO NOT CHANGE THIS LINE, fix problems instead.
    app.warningiserror = False


intersphinx_mapping = dict(
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    xarray=("http://xarray.pydata.org/en/stable/", None),
)
qualname_overrides = {
    "anndata._core.anndata.AnnData": "anndata.AnnData",
}
