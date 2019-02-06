from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='sbio_project',
    version='1.0',
    description='project to compute geometric feature of RepeatDB protein',
    # url
    author='Federico Baldo',
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Programming Language :: Python :: 2.7'
    ],
    author_email='federico.baldo.1@studenti.unipd.it',
    install_requires=["biopython", "numpy", "pymol", "pyinquirer"],
    python_requires='>=2.7',
    entry_points={  
        'console_scripts': [
            'main=bio_project.cli:main',
        ],
    }
)