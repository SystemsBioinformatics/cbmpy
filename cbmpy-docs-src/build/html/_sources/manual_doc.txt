.. _Introducing-CBMPy:

Introduction
============

PySCeS CBMPy is a new platform for constraint based modelling and analysis. It has been designed using principles
developed in the PySCeS simulation software project: usability, flexibility and accessibility. Its architecture 
is both extensible and flexible using data structures that are intuitive to the biologist (metabolites, reactions, compartments)
while transparently translating these into the underlying mathematical structures used in advanced analysis (LP's, MILP's).

PySCeS CBMPy implements popular analyses such as FBA, FVA, element/charge balancing, network analysis and model editing as
well as advanced methods developed specifically for the ecosystem modelling: minimal distance methods, flux minimization and input selection.

To cater for a diverse range of modelling needs PySCeS CBMPy supports user interaction via:

- interactive console, scripting for advanced use or as a library for software development
- GUI, for quick access to a visual representation of the model, analysis methods and annotation tools
- SOAP based web services: using the Mariner framework much high level functionality is exposed for integration into web tools

For more information on the development and use of PySCeS CBMPy visit the website (http:cbmpy.sourceforge.net) for up to date information and 
feel free to contact the development team (bgoli@users.sourceforge.net).
