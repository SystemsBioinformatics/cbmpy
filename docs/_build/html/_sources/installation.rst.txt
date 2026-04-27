Installation Guide
==================

Installing CBMPy
----------------

CBMPy supports Python 3.8 through 3.12. You can install it using conda or pip.

Prerequisites
-------------

**System Requirements:**

- Python 3.8 - 3.12
- One of the following solvers:
  - IBM CPLEX (commercial)
  - GLPK (open source, via `swiglpk`)
- numpy >= 1.17.0
- python-libsbml >= 5.18.0 (for SBML support)
- scipy (for numerical operations)

Using Conda (Recommended)
-------------------------

1. Create a conda environment:

.. code-block:: bash

   conda create -n cbmpy3 python=3.10

2. Activate the environment:

.. code-block:: bash

   conda activate cbmpy3

3. Install CBMPy:

.. code-block:: bash

   conda install cbmpy

Alternatively, create an environment with all dependencies:

.. code-block:: bash

   conda env create -f environment.yml

Using pip
---------

1. Create a virtual environment:

.. code-block:: bash

   python -m venv cbmpy3
   source cbmpy3/bin/activate  # On Linux/Mac
   cbmpy3\Scripts\activate.bat  # On Windows

2. Install CBMPy:

.. code-block:: bash

   pip install -e ".[dev]"

Or from the requirements file:

.. code-block:: bash

   pip install -r requirements.txt

Verification
------------

Verify the installation by running:

.. code-block:: python

   import cbmpy
   cmod = cbmpy.readSBML3FBC('cbmpy_test_core')
   cbmpy.doFBA(cmod)

Troubleshooting
---------------

**Error: "No linear solver present"**

You need to install either CPLEX or GLPK. On Linux:

.. code-block:: bash

   # Install GLPK via swiglpk
   pip install swiglpk

**Error: "python-libsbml not found"**

Install python-libsbml:

.. code-block:: bash

   pip install python-libsbml

**On Windows with CPLEX:**

1. Download the CPLEX Python package from IBM
2. Add the CPLEX directory to your PYTHONPATH
3. Install CBMPy

Update CBMPy
------------

After installation, you can update CBMPy using:

.. code-block:: bash

   # Conda
   conda update cbmpy

   # pip
   pip install --upgrade cbmpy

Or update in-place:

.. code-block:: bash

   pip install -e ".[dev]"

Dependencies
------------

Core dependencies:
- numpy >= 1.17.0
- packaging
- pyparsing
- python-libsbml >= 5.18.0
- lxml
- xlrd
- xlwt
- swiglpk
- scipy
- XlsxWriter

Optional dependencies:
- sympy >= 1.7.1 (for symbolic matrix support)

Development dependencies:
- pytest >= 7.0
- pytest-cov
- docformatter >= 1.7.7
- sphinx >= 7.0
- sphinx_rtd_theme
- myst-parser
- tox
