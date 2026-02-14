# Configuration file for the Sphinx documentation builder.
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

project = "NCHASH"
copyright = "2024, NCHASH Authors"
author = "NCHASH Authors"

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

# Include docs directory
import os
import sys
sys.path.insert(0, os.path.abspath(".."))
