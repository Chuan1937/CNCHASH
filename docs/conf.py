# Configuration file for the Sphinx documentation builder.
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Include docs directory
sys.path.insert(0, os.path.abspath(".."))

project = "CNCHASH"
copyright = "2026, He XingChen"
author = "He XingChen"

# Extensions
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# MyST parser settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

# Source files
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# HTML theme
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Last updated format
html_last_updated_fmt = "%Y-%m-%d"

# Extra context for templates
html_context = {
    "display_author": True,
}
