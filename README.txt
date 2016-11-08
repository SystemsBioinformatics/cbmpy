CBMPy
=====

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
Last edit: $Author: bgoli $ ($Id: CBConfig.py 144 2014-01-16 22:27:55Z bgoli $)

Installation
============

PySCeS CBMPy (http://cbmpy.sourceforge.net) is a new platform for constraint
based modelling and analysis. It has been designed using principles developed
in the PySCeS simulation software project: usability, flexibility and accessibility. 
Its architecture is both extensible and flexible using data structures that are intuitive 
to  the biologist (metabolites, reactions, compartments) while transparently translating 
these into the underlying mathematical structures used in advanced analysis (LP's, MILP's).

PySCeS CBMPy implements popular analyses such as FBA, FVA, element/charge
balancing, network analysis and model editing as well as advanced methods
developed specifically for the ecosystem modelling: minimal distance methods,
flux minimization and input selection. To cater for a diverse range of modelling
needs PySCeS CBMPy supports user interaction via:

- interactive console, scripting for advanced use or as a library for software development
- GUI, for quick access to a visual representation of the model, analysis methods and annotation tools
- SOAP based web services: using the Mariner framework much high level functionality is exposed for integration into web tools

For more information on the development and use of PySCeS CBMPy feel free to contact me:

PySCeS-CBMPy has been tested on Windows 7 and 8.1, Mac OSX and Ubuntu Linux 12.04, 14.04, 16.04. 
It is compatible with both Python 2.7+ and includes experimental support for Python 3.4+ It is highly recommend to use
Python 2.7 as not all Python package dependencies (extended functionality) are available for Python 3.

PySCeS CBMPy is now accessible as a Python module **cbmpy** in place of the the previously used **pyscescbm** which is no longer supported. CBMPy includes support for  reading/writing models in SBML3 FBC versions 1 and 2 as well as COBRA dialect, Excel spreadsheets and Python.

To use follow the installation instructions given below and try the following in a Python shell::

import cbmpy
cmod = cbmpy.readSBML3FBC('cbmpy_test_core')
cbmpy.doFBA(cmod)

New Ipython notebook tutorials are available. Happy modelling!

New cross-platform conda support
--------------------------------

CBMPy is now available as a crosss platform Anaconda package available for Win32/64, Linux32/64 and OSX. This is the first
major release using the conda platform and all feedback is welcomed. If you have Anaconda installed (tested with 4.2.0)
simply try the following::

 conda install -c bgoli -c SBMLTeam cbmpy
 
which should install CBMPy, libSBML and all the required (non-solver) dependencies. If you do not use Anaconda or would
like to try install things in a semi-automated way carry on reading.

The following installation instructions are for Ubuntu 16.04 but should be adaptable to any
Linux package managment system, OSX, Debian, etc. Except for GLPK (4.47) and SymPy (0.7.4 or newer)
no specific library version is required. For more detailed installation instructions and Windows
please see the online documentation http://cbmpy.sourceforge.net/reference/install_doc.html

Auto-dependency configuration
-----------------------------

I am in the process of creating automated dependency checking and building tools for CBMPy. 

Ubuntu support is complete with cross platform conda support now available for testing (https://github.com/bgoli/cbmpy-build.git)::

To start with download the install script that will run the system test::

 curl --remote-name https://raw.githubusercontent.com/bgoli/cbmpy-build/master/ubuntu/1_install_cbmpy_dependencies.py 

This script will analyse your Python distribution for missing dependencies and generate a shell script that will auto-install them.

For an installation using the native Python supplied with Ubuntu run::

 python 1_install_cbmpy_dependencies.py UBUNTU
 
For a local user installation of Anaconda try the following. This is not strictly necessary if only CBMPy is going to be installed (see above)::

 python 1_install_cbmpy_dependencies.py CONDA

in both cases a shell script `metatoolkit_install_dependencies.sh` is generated and printed to the screen. Run this script 
or install the relevant packages by hand.

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

 conda install -c bgoli cbmpy

or download the source  and run::

 python setup.py build sdist
 sudo python setup.py install

Python3 (experimental)
----------------------

Not all dependencies are available for Python3::

sudo apt-get install python3-dev python3-numpy python3-scipy python3-matplotlib  python3-pip
sudo apt-get install python3-xlrd python3-h5py

# need to find out what is going on with Python3 and xlwt suds
# easy_install3 sympy ???
# wxPython and PyQt4 not in Ubuntu P3 builds yet

sudo apt-get install ipython3 ipython3-notebook

sudo apt-get install libxml2 libxml2-dev
sudo apt-get install zlib1g zlib1g-dev
sudo apt-get install bzip2 libbz2-dev

sudo pip3 install python-libsbml-experimental

sudo apt-get install python-sip python-sip-dev build-essential

More information in the docs/ directory.

(c) Brett G. Olivier, Amsterdam, 2014-2016