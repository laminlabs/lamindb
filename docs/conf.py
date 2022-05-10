import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE)]
from lamin_sphinx import *  # noqa
import lamindb  # noqa

project = "lamindb"
html_title = f"{project} | Lamin Labs"
release = lamindb.__version__
html_context["github_repo"] = "lamindb"  # noqa
html_sidebars = {
    "*": ["sidebar-nav-bs"],
    "**/*": ["sidebar-nav-bs"],
}
