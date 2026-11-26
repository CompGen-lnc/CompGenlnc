from dataclasses import fields

from docutils import nodes
from sphinx import addnodes

project = "CompGenlnc"
copyright = "2026, Cindy Calderón Arce, Gabriel Calvo Montero"
author = "Cindy Calderón Arce, Gabriel Calvo Montero"
release = "0.1.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "numpydoc",
    "rinoh.frontend.sphinx",
]

templates_path = ["_templates"]
exclude_patterns = []
default_role = 'literal'

# Autodoc Settings
autodoc_class_signature = "separated"
autodoc_default_options = {
    "inherited-members": False,
    "exclude-members": "__*",
    "member-order": "bysource",
    "no-index": True,
    "signature": False,
}
autodoc_preserve_defaults = True

# HTML Settings
html_theme = "furo"
html_static_path = ["_static"]

rinoh_documents = [
    dict(
        doc='index',
        target='compgenlnc_manual',
        title=project,
        subtitle=release,
        author=author,
        template='rinohtype_template.rtt',
    )
]

module_alias = {
    ".consts.paths": ".paths",
    ".consts.fields": "",
    ".consts.params": "",
    ".consts.regex": "",
    ".dtypes": "",
}


def change_module_alias(app, doctree, docname):
    for node in doctree.traverse(addnodes.desc_addname):
        text = node.astext()
        for module in module_alias:
            if module in text:
                new_text = text.replace(module, module_alias[module])
                node.children = [nodes.Text(new_text)]


def setup(app):
    app.connect('doctree-resolved', change_module_alias)
