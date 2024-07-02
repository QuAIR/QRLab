
import os
import sys

sys.path.insert(0, os.path.join("..", "..", ".."))
matlab_src_dir = os.path.abspath(os.path.join("..", "..", ".."))
extensions = [
    "sphinx.ext.autodoc",
    "sphinxcontrib.matlab",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx_immaterial",
]
primary_domain = "mat"
# project = "QRLab"
master_doc = "index"
source_suffix = ".rst"
nitpicky = True


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "sphinx_rtd_theme"
html_theme = "sphinx_immaterial"
html_title = "QRLab"
html_short_title = "QRLab"
build_dir = "api"
# html_theme_options = {
#     'navigation_depth': 1,
# }
html_theme_options = {
    "base_url": "https://quair.github.io/QRLab/",
    "repo_url": "https://github.com/QuAIR/QRLab",
    "repo_name": "QRLab",
    # 'google_analytics_account': 'UA-XXXXX',
    "html_minify": True,
    "css_minify": True,
    "nav_title": "QRLab API Documentation",
    # 'logo_icon': '&#xe869',
    # 'globaltoc_depth': 2,
    "palette": { "primary": "orange" }
}
html_favicon = '../favicon.svg'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

master_doc = "index"

# Autodoc
napoleon_numpy_docstring = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_warningiserror = False
autodoc_inherit_docstrings = False
autodoc_docstring_signature = False
autodoc_typehints_description_target = "documented"
autodoc_typehints_format = "short"
