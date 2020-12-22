import sys

extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'PeleC'
copyright = u'AMReX Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory and the Alliance for Sustainable Energy, LLC., through National Renewable Energy Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.'
author = u'J.B. Bell, M.S. Day, E. Motheau, D. Graves, M. Henry de Frahan, R.W. Grout, J. Rood'
version = u'0.1'
release = u'0.1'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
htmlhelp_basename = 'PeleCdoc'
latex_elements = { }
latex_documents = [
    (master_doc, 'PeleC.tex', u'PeleC Documentation',
     author, 'manual'),
]
texinfo_documents = [
    (master_doc, 'PeleC', u'PeleC Documentation',
     author, 'PeleC', 'One line description of project.',
     'Miscellaneous'),
]
