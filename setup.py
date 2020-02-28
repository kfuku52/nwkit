from setuptools import setup, find_packages

setup(
        name             = 'nwkit',
        version          = "0.2",
        description      = 'Tools for processing newick trees',
        license          = "BSD 3-clause License",
        author           = "Kenji Fukushima",
        author_email     = 'kfuku52@gmail.com',
        url              = 'https://github.com/kfuku52/nwkit.git',
        keywords         = '',
        packages         = find_packages(),
        install_requires = ['ete3',],
        scripts          = ['nwkit/nwkit',],
)