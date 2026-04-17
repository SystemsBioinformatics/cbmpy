CBMPy Documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   api/modules
   api/core
   fluxmodules

.. note::

   CBMPy is a platform for constraint-based modelling and analysis of biological systems.
   It implements Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), and other
   constraint-based methods using linear and mixed-integer programming.

.. toctree::
   :maxdepth: 2
   :hidden:

   release-notes

Installation Guide
------------------

See :doc:`installation` for installation instructions.

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

For more examples, see the :doc:`installation` guide.

API Documentation
-----------------

.. toctree::
   :maxdepth: 2

   api/modules
   api/core
   fluxmodules

Release Notes
-------------

.. toctree::
   :maxdepth: 1

   release-notes

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
