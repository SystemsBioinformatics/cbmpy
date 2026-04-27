CBMPy Documentation
===================

CBMPy is a platform for constraint-based modelling and analysis of biological systems.
It implements Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), and other
constraint-based methods using linear and mixed-integer programming.

Quick Start
-----------

After installing CBMPy:

.. code-block:: python

   import cbmpy

   # Read a model from SBML
   cmod = cbmpy.readSBML3FBC('example_model.xml')

   # Perform Flux Balance Analysis
   result = cbmpy.doFBA(cmod)

   # Print the optimal objective value
   print(f"Optimal flux: {result[0]}")

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   api/modules
   api/core
   api/solvers
   api/io
   api/networking
   api/visualization
   release-notes

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
