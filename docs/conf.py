"""Sphinx configuration."""
project = "Ont_Chopper"
author = "Ting-You Wang"
copyright = "2022, Ting-You Wang"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "myst_parser",
]
autodoc_typehints = "description"
html_theme = "furo"
