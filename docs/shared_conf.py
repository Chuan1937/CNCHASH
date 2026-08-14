# Shared Sphinx configuration for CNCHASH docs.
# docs/en/conf.py and docs/zh_CN/conf.py import from here and set the
# language. Follows the layout used by efwi3D.

import os
import subprocess
from datetime import datetime
from pathlib import Path


project = "CNCHASH"
author = "He XingChen"
copyright = f"{datetime.now().year}, {author}"

repository_root = Path(__file__).resolve().parents[1]

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

templates_path = ["../_templates"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["../_static"]
html_last_updated_fmt = "%Y-%m-%d"


def _git_last_commit(path: Path) -> tuple[str, str]:
    """(commit date, short hash) of the last commit touching a file."""
    result = subprocess.run(
        ["git", "log", "-1", "--format=%ci%x00%h", "--", str(path)],
        cwd=repository_root,
        capture_output=True,
        text=True,
        timeout=5,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return "", ""
    date, commit = result.stdout.strip().split("\x00", maxsplit=1)
    return date, commit


def setup(app):
    """Show the last git commit date of each page in the footer."""
    def html_page_context_handler(app, pagename, templatename, context, doctree):
        suffixes = list(app.config.source_suffix.keys()) if not isinstance(
            app.config.source_suffix, list
        ) else app.config.source_suffix
        for suffix in suffixes:
            source_file = Path(app.srcdir) / (pagename + suffix)
            if source_file.exists():
                date, commit = _git_last_commit(source_file)
                if date:
                    context["last_updated"] = date[:10]
                if commit:
                    context["commit_hash"] = commit
                break
        context["author"] = app.config.author

    app.connect("html-page-context", html_page_context_handler)
