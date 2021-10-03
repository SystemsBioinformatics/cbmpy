"""
CBMPy: setup.py
===============
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2010-2022 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: setup.py 698 2019-07-30 14:03:26Z bgoli $)

"""

import os
import time

local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))

try:
    from setuptools import setup

    install_requires_src = ['numpy', 'packaging', 'nose']
    extras_require_src = {
        'sympy': ['sympy'],
        'glpk': ['swiglpk',],
        'sbml': ['python_libsbml', 'lxml',],
        'all': ['sympy', 'swiglpk', 'python_libsbml', 'lxml',],
    }
    tests_require_src = ['nose']
except:
    from distutils.core import setup

    install_requires_src = []
    extras_require_src = {}
    tests_require_src = []

mydata_files = []

# release get version
vmaj = vmin = vmic = ''
with open(os.path.join(local_path, 'cbmpy', 'CBConfig.py')) as F:
    for l in F:
        if l.startswith('__VERSION_MAJOR__'):
            vmaj = l.split('=')[1].strip()
        elif l.startswith('__VERSION_MINOR__'):
            vmin = l.split('=')[1].strip()
        elif l.startswith('__VERSION_MICRO__'):
            vmic = l.split('=')[1].strip()

mypackages = ['cbmpy', 'cbmpy.fluxmodules']

readme = """ CBMPy
 =====

 PySCeS CBMPy (http://cbmpy.sourceforge.net) is a new platform for constraint
 based modelling and analysis. It has been designed using principles developed
 in the PySCeS simulation software project: usability, flexibility and accessibility. Its architecture is both extensible and flexible using data structures that are intuitive to  the biologist (metabolites, reactions, compartments) while transparently translating these into the underlying mathematical structures used in advanced analysis (LP's, MILP's).

 PySCeS CBMPy implements popular analyses such as FBA, FVA, element/charge
 balancing, network analysis and model editing as well as advanced methods
 developed specifically for the ecosystem modelling: minimal distance methods,
 flux minimization and input selection. To cater for a diverse range of modelling
 needs PySCeS CBMPy supports user interaction via:

 - interactive console, scripting for advanced use or as a library for software development
 - GUI, for quick access to a visual representation of the model, analysis methods and annotation tools
 - SOAP based web services: using the Mariner framework much high level functionality is exposed for integration into web tools

 For more information on the development and use of PySCeS CBMPy feel free to contact me:

 PySCeS-CBMPy has been tested on Windows 7 and 8.1, Mac OSX and Ubuntu Linux 12.04, 14.04, 16.04. It is compatible with Python 3.6+ 
 
  To use follow the installation instructions given below and try the following in a Python shell::

 import cbmpy
 cmod = cbmpy.readSBML3FBC('cbmpy_test_core')
 cbmpy.doFBA(cmod)

"""

setup(
    package_dir={'cbmpy': 'cbmpy'},
    packages=mypackages,
    data_files=mydata_files,
    long_description=readme,
    name="cbmpy",
    summary="CBMPy: Constraint Based Modelling in Python",
    version='{}.{}.{}'.format(vmaj, vmin, vmic),
    maintainer='Brett G. Olivier',
    author='Brett G. Olivier',
    author_email='b.g.olivier@vu.nl',
    maintainer_email='b.g.olivier@vu.nl',
    url="https://github.com/SystemsBioinformatics/cbmpy",
    download_url="https://github.com/SystemsBioinformatics/cbmpy/releases",
    license="GNU General Public License (GPL)",
    keywords="computational systems biology, modelling, simulation, genome scale models, sbml, constraint-based modelling, fbc, linear programming, groups, standard",
    # zip_safe = False,
    install_requires=install_requires_src,
    extras_require=extras_require_src,
    tests_require=tests_require_src,    
    platforms=["Windows", "Linux", "Mac"],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    description=""" CBMPy (https://github.com/SystemsBioinformatics/cbmpy) is a platform for constraint based modelling and analysis. CBMPy implements popular analyses such as FBA, FVA, element/charge balancing, network analysis and model editing as well as advanced methods developed specifically for the ecosystem modelling.""",
)
try:
    os.chdir(os.path.join(local_path, 'cbmpy', 'fluxmodules'))
    print(os.getcwd())
    import cbmpy

    print(cbmpy.__version__)
    os.chdir(local_path)
except ImportError:
    pass

