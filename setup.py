from setuptools import setup, find_packages
import os
import re
import ast

with open(os.path.join('nwkit', '__init__.py')) as f:
        match = re.search(r'__version__\s+=\s+(.*)', f.read())
version = str(ast.literal_eval(match.group(1)))

setup(
        name             = 'nwkit',
        version          = version,
        description      = 'Tools for processing newick trees',
        license          = "BSD 3-clause License",
        author           = "Kenji Fukushima",
        author_email     = 'kfuku52@gmail.com',
        url              = 'https://github.com/kfuku52/nwkit.git',
        keywords         = 'phylogenetics',
        packages         = find_packages(),
        install_requires = ['ete3','biopython','pandas','requests','numpy',],
        scripts          = ['nwkit/nwkit',],
        include_package_data = True,
        package_data     = {
                '':['data_tree/*.nwk',],
        }
)