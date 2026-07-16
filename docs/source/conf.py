# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# -- Project information
project = 'AstroGeoFit'
copyright = '2025, CNRS -- Contact: astrogeofit@astrogeo.eu'
author = 'Jacques Laskar et al. (CNRS)'
release = '1.0'

# -- General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_togglebutton',
]

templates_path = ['_templates']
html_static_path = ['_static']
html_theme = 'furo'

html_css_files = ['custom.css']

# Your custom footer content
html_context = {
    "custom_footer": "2025, CNRS<br>Contact: <a href='mailto:astrogeofit@astrogeo.eu'>astrogeofit@astrogeo.eu</a>"
}

