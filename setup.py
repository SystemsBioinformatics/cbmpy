"""
CBMPy: setup.py
===============
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2010-2016 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: setup.py 572 2017-03-26 13:00:08Z bgoli $)

"""

import os
local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))

USE_SETUPTOOLS = False

if USE_SETUPTOOLS:
    # setuptools
    from setuptools import setup
else:
    # distutils
    from distutils.core import setup

# new more elegant hack implemented to get rid of evil hack and allow pip install
"""
# I suspect this is a distutils specific evil hack to get the data into the install directory
from distutils.command.install import INSTALL_SCHEMES
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']
"""
mydata_files = []
# new hack should get rid of the need for this
#mydata_files.append((os.path.join('cbmpy', 'models'),\
#[os.path.join(local_path,'src','models',examplefile) for examplefile in os.listdir(os.path.join(local_path,'src','models'))\
#if examplefile.endswith('.xml')]))

# new hack should get rid of the need for this
#mydata_files.append((os.path.join('cbmpy', 'nosetests'),\
#[os.path.join(local_path,'src','nosetests',examplefile) for examplefile in os.listdir(os.path.join(local_path,'src','nosetests'))\
#if examplefile.endswith('.xml') or examplefile.endswith('.json')]))

# release
try:
    STATUS = 'r'+'$Rev: 572 $'.replace('Rev: ','').replace('$','').strip()
except:
    STATUS = 'beta'

mypackages= ['cbmpy', 'cbmpy.fluxmodules', 'cbmpy.models', 'cbmpy.nosetests', 'cbmpy.solver']

setup(
 package_dir = {'cbmpy' : 'src'},
 packages = mypackages,
 data_files = mydata_files,
 name = "cbmpy",
 summary = "CBMPy: Constraint Based Modelling in Python",
 #version = "0.7.12.{}".format(STATUS),
 version = "0.7.12",
 maintainer = 'Brett G. Olivier',
 author = 'Brett G. Olivier',
 maintainer_email = 'bgoli@users.sourceforge.net',
    url = "http://cbmpy.sourceforge.net",
    download_url = "http://sourceforge.net/projects/cbmpy/files/release/0.7/",
    license = "GNU General Public License (GPL)",
    keywords = "computational systems biology, modelling, simulation, genome scale models, sbml, constraint-based modelling, fbc, linear programming, groups, standard" ,
    #zip_safe = False,
    requires = ['numpy', 'sympy', 'libsbml', 'nose'],
    platforms = ["Windows", "Linux", "Mac"],
    classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry'
     ],
 description =\
 """ PySCeS-CBMPy
 ============

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

 PySCeS-CBMPy has been tested on Windows 7 and 8.1, Mac OSX and Ubuntu Linux 12.04, 14.04, 16.04. It is compatible with both Python 2.7+ and includes experimental support for Python 3.4+ It is highly recommend to use
 Python 2.7 as not all Python package dependencies (extended functionality) are available for Python 3.

 PySCeS CBMPy is now accessible as a Python module **cbmpy** in place of the the previously used **pyscescbm** which is no longer supported. CBMPy includes support for  reading/writing models in SBML3 FBC versions 1 and 2 as well as COBRA dialect, Excel spreadsheets and Python.

 To use follow the installation instructions given below and try the following in a Python shell::

 import cbmpy
 cmod = cbmpy.readSBML3FBC('cbmpy_test_core')
 cbmpy.doFBA(cmod)

 New Ipython notebook tutorials are available. Happy modelling!

 The following installation instructions are for Ubuntu 16.04 but should be adaptable to any
 Linux package managment system, OSX, Debian, etc. Except for GLPK (4.47) and SymPy (0.7.4 or newer)
 no specific library version is required. For more detailed installation instructions and Windows
 please see the online documentation http://cbmpy.sourceforge.net/reference/install_doc.html

 New! auto-dependency configuration
 ----------------------------------

 I am in the process of creating automated dependency checking and building tools for CBMPy. These can be found at::

  https://github.com/bgoli/cbmpy-build

 Ubuntu support is almost complete with Windows/Conda support in development, grab form GitHub::

  https://github.com/bgoli/cbmpy-build.git

 Manual dependency configuration is provided below. For Windows users most of these utilities are included in
 Python distributions like Anaconda (recomended)

 Python2
 -------

 First we create a scientific Python workbench::

 sudo apt-get install python-dev python-numpy python-scipy python-matplotlib  python-pip
 sudo apt-get install python-sympy python-suds python-xlrd python-xlwt python-h5py
 sudo apt-get install python-wxgtk2.8
 sudo apt-get install ipython ipython-notebook

 libSBML
 ~~~~~~~

 Installing libSBML is now easy using Pip::

 sudo apt-get install libxml2 libxml2-dev
 sudo apt-get install zlib1g zlib1g-dev
 sudo apt-get install bzip2 libbz2-dev

 sudo pip install --update python-libsbml

 Extended functionality
 ~~~~~~~~~~~~~~~~~~~~~~

 sudo pip install biopython docx

 Windows
 ~~~~~~~

 Use easy_install, pip or your package manager (e.g. conda) to install the following packages::

  numpy scipy matplotlib sympy xlrd xlwt
  biopython docx suds wxPython

 pip install --update python-libsbml

 glpk/python-glpk
 ~~~~~~~~~~~~~~~~

 CBMPy requires a linear solver for numerical analysis, the open source (glpk) solver can be automatically built and installed as follows (requires git to be installed and accessible):

 Download the install script that will install GLPK/PyGLPK for CBMPy on Ubuntu 14.04 or newer::

  curl --remote-name https://raw.githubusercontent.com/bgoli/cbmpy-glpk/master/install_glpk.sh

 Make executable::

  chmod 744 install_glpk.sh

 and run::

  ./install_glpk.sh

 Note this script is designed to be used for building containers and will remove any installed version of GLPK and build and install the correct version needed for PyGLPK.

 No warranty of any kind assumed or otherwise, use at own risk!

 CBMPy
 ~~~~~

 Finally, install CBMPy::

  sudo easy_install cbmpy

 or

  sudo pip install cbmpy

 or try the new experimental CONDA support::

  conda install -c bgoli -c sbmlteam cbmpy

 or download the source  and run::

  python setup.py build sdist
  sudo python setup.py install

""")
try:
    import cbmpy
except ImportError:
    pass
