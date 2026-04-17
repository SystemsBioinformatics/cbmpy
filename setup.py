"""
CBMPy: setup.py
===============
PySCeS Constraint Based Modelling (https://systemsbioinformatics.github.io/cbmpy/)

This setup.py reads from pyproject.toml for pyproject.toml-based builds.
This file is kept for backward compatibility with legacy setup.py tools.
"""

from setuptools import setup

# Read version from CBConfig
import os
import sys

local_path = os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), 'cbmpy', 'CBConfig.py')))

try:
    from setuptools import setup

    # Read version from CBConfig
    vmaj = vmin = vmic = ''
    with open(os.path.join(local_path, 'cbmpy', 'CBConfig.py')) as F:
        for l in F:
            if l.startswith('__VERSION_MAJOR__'):
                vmaj = l.split('=')[1].strip()
            elif l.startswith('__VERSION_MINOR__'):
                vmin = l.split('=')[1].strip()
            elif l.startswith('__VERSION_MICRO__'):
                vmic = l.split('=')[1].strip()

    install_requires_src = ['numpy', 'packaging', 'pyparsing', 'python_libsbml', 'lxml', 'xlrd', 'xlwt', 'swiglpk', 'scipy', 'XlsxWriter']
    extras_require_src = {
        'all': ['sympy', 'numpy', 'packaging', 'pyparsing', 'python_libsbml', 'lxml', 'xlrd', 'xlwt', 'swiglpk', 'scipy', 'nose_py3', 'XlsxWriter'],
    }
    tests_require_src = ['numpy', 'packaging', 'pyparsing', 'python_libsbml', 'lxml', 'numpy', 'nose_py3']

    mydata_files = []
    mypackages = ['cbmpy', 'cbmpy.fluxmodules']

    description = """CBMPy: https://systemsbioinformatics.github.io/cbmpy/ is a platform for constraint based modelling. It implements popular FBA analyses, model editing and supports the latest SBML3FBC standards."""

    setup(
        package_dir={'cbmpy': 'cbmpy'},
        packages=mypackages,
        data_files=mydata_files,
        long_description=description,
        name="cbmpy",
        summary="CBMPy: Constraint Based Modelling in Python",
        version='{}.{}.{}'.format(vmaj, vmin, vmic),
        maintainer='Brett G. Olivier',
        author='Brett G. Olivier',
        author_email='b.g.olivier@vu.nl',
        maintainer_email='b.g.olivier@vu.nl',
        url="https://systemsbioinformatics.github.io/cbmpy/",
        download_url="https://github.com/SystemsBioinformatics/cbmpy/releases",
        license="GNU General Public License (GPL)",
        keywords="computational systems biology, modelling, simulation, genome scale models, sbml, constraint-based modelling, fbc, linear programming",
        install_requires=install_requires_src,
        extras_require=extras_require_src,
        tests_require=tests_require_src,
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Environment :: MacOS X',
            'Environment :: Win32 (MS Windows)',
            'Environment :: X11 Applications',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Software Development :: Libraries :: Python Modules',
        ],
        description=description,
    )
except Exception as e:
    print(f"Error reading version from CBConfig: {e}")
    from distutils.core import setup
    setup()
