from mutalyzer_crossmapper import _get_metadata


author = _get_metadata('Author')
copyright = _get_metadata('Author')
project = _get_metadata('Name')
release = _get_metadata('Version')

autoclass_content = 'both'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx'
    ]
master_doc = 'index'
