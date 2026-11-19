# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "CompGenlnc"
copyright = "2026, Cindy Calderón Arce, Gabriel Calvo Montero"
author = "Cindy Calderón Arce, Gabriel Calvo Montero"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = []

# Autodoc Settings
autodoc_class_signature = "separated"
autodoc_default_options = {
    "inherited-members": False,
    "exclude-members": "__*",
    "member-order": "bysource",
    "no-index-entry": True,
    "noindex": True,
    "signature": False,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
