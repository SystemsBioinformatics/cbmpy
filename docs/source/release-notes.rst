Release Notes
-=-=--=-=----

0.8.9 (unreleased)
--=-=---------

Docstring improvements:

- Converted all module docstrings to numpy format across the entire codebase
  - CBXML.py: 68 docstrings converted
  - CBModel.py: 206 docstrings converted
  - CBCPLEX.py: 25 docstrings converted
  - CBGLPK.py: 16 docstrings converted
  - CBNetDB.py: 21 docstrings converted
  - CBRead.py: 6 docstrings converted
  - CBWrite.py: 31 docstrings converted
- All docstrings now use numpy format compatible with Sphinx Napoleon extension
- Sphinx documentation updated with all documented modules

Improvements:

- Updated Sphinx docs with new module API pages
  - api/solvers.rst - CBSolver, CBSolverX, CBCPLEX, CBGLPK, CBGLPKOLD
  - api/io.rst - CBRead, CBReadtxt, CBWrite, CBXML
  - api/networking.rst - CBNetDB, miriamids
  - api/visualization.rst - CBPlot, GUI frontends, multi-processing, utilities
- Added missing core classes to docs (Species, Objective, GeneProteinAssociation, UserConstraint, Group)
- Fixed docstring style reference in CLAUDE.md (numpydoc -> numpy)

0.8.8 (2023)
-==-=-=---

- Support for Python 3.8 through 3.12
- SBML FBCv3 support
- Flux variability analysis improvements

0.8.4 (2022)
-==-=-=--

- Python 3.x compatibility update
- SBML3 FBC standard support improvements

0.8.2 (2021)
-==-=-=--

- Bug fixes for SymPy compatibility
- Improved requirements handling
