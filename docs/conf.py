import os
import subprocess
import sys

subprocess.call('pip install ..', shell=True)

from crossmapper import _get_metadata


author = _get_metadata('Author')
copyright = _get_metadata('Author')
project = _get_metadata('Name')
release = _get_metadata('Version')

autoclass_content = 'both'
extensions = ['sphinx.ext.autodoc']
master_doc = 'index'
