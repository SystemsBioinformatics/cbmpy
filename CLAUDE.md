# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# CBMPy - Constraint Based Modelling in Python

CBMPy is a Python platform for constraint-based modelling and analysis of biological systems. It implements Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), and other constraint-based methods using linear and mixed-integer programming.

## Installation

**Use conda environments only** - legacy pip installation is deprecated.

```bash
conda env create -f=environment.yml
conda activate cbmpy3
```

After installation, verify:

```python
import cbmpy
cmod = cbmpy.readSBML3FBC('cbmpy_test_core')
cbmpy.doFBA(cmod)
```

## Development Environment

```bash
# Activate the environment
conda activate cbmpy3

# Install dev dependencies
pip install -e ".[dev]"

# Run the test suite
pytest
```

## Common Development Tasks

### Running tests
```bash
# Install pytest
pip install pytest pytest-cov

# Run full test suite
pytest

# Run a specific test file
pytest cbmpy/tests/test_specific.py

# Run a specific test function
pytest cbmpy/tests/test_specific.py::test_something

# Run tests with verbose output and coverage
pytest -v --cov=cbmpy

# Run tests with verbose output and coverage report
pytest -v --cov=cbmpy --cov-report=html

# Run tests excluding slow tests
pytest -m "not slow"

# Run only integration tests
pytest -m integration

### Building and packaging
```bash
# Build source distribution
python setup.py sdist

# Build wheel
python setup.py bdist_wheel

# Clean and rebuild
rm -rf build dist *.egg-info
python setup.py sdist bdist_wheel
```

### Code style and linting
The project uses pylint with custom disable comments for:
- C0103 (invalid variable name)
- C0301 (line too long)
- E1101 (module has no member)

Code is formatted with standard Python conventions (no AI auto-formatting).

### Docstring formatting

All docstrings must follow the **numpydoc** format to be compatible with Sphinx documentation generation.

**Docstring structure:**

```python
def function_name(param1, param2):
    """
    Short description of the function.

    Extended description (if needed).

    Parameters
    ----------
    param1 : type
        Description of param1.
    param2 : type
        Description of param2.

    Returns
    -------
    return_type
        Description of return value.

    Raises
    ------
    ExceptionType
        Description of when and why this exception is raised.

    Notes
    -----
    Any important notes or warnings.

    Examples
    --------
    >>> example_code_here()
    """
```

**Formatting commands:**

```bash
# Format all docstrings to numpydoc style
docformatter --in-place --recursive --style=numpy cbmpy/

# Check docstring formatting (should have no output)
docformatter --diff --recursive cbmpy/

# Format specific file
docformatter --in-place cbmpy/module_name.py
```

### Building and packaging
```bash
# Build source distribution
python setup.py sdist

# Build wheel
python setup.py bdist_wheel

# Clean and rebuild
rm -rf build dist *.egg-info
python setup.py sdist bdist_wheel
```

### Code style and linting
The project uses pylint with custom disable comments for:
- C0103 (invalid variable name)
- C0301 (line too long)
- E1101 (module has no member)

Code is formatted with standard Python conventions (no AI auto-formatting).

## Project Architecture

### Main Package: `cbmpy/`

**Core modules:**

- **CBModel.py** - Core model class representing metabolic networks
  - Contains `Fbase` base class for all model objects
  - Manages metabolites, reactions, and compartments
  - Provides SBML read/write functionality

- **CBSolver.py** - Solver implementations (CPLEX, GLPK)
  - Implements LP/MILP solvers for FBA, FVA
  - Main entry points: `analyzeModel`, `FluxVariabilityAnalysis`

- **CBRead.py / CBWrite.py** - File I/O
  - `readSBML3FBC`, `readSBML2FBA`, `readCOBRASBML`
  - `writeSBML3FBC`, `writeSBML3FBCV3` (FBC v3 support)

- **CBCommon.py** - Common utilities
  - Matrix operations (`StructMatrixLP`)
  - MIRIAM annotation parsing
  - ID checking and fixing

- **CBSolverX.py** - Extended solver (QP support)
  - Handles quadratic objectives
  - Additional constraint types

**Supporting modules:**

- **CBConfig.py** - Configuration and version info
- **CBXML.py** - SBML handling (uses python-libsbml)
- **CBNetDB.py** - Network analysis
- **CBPlot.py** - Visualization
- **CBTools.py** - Utility functions
- **PyscesSED.py** - SED-ML support

### Flux Modules: `cbmpy/fluxmodules/`

Advanced analyses for module decomposition:

- **fluxmodules.py** - Main flux module interface
  - `computeModules()` - Compute flux modules
  - `computeModulesMatroid()` - Matroid-based computation
  - Requires prior FVA to identify variable reactions

- **matroid.py** - Matroid theory implementation
- **sparserationals.py** - Sparse rational number arithmetic
- **enumeration.py** - Module enumeration algorithms

- **Test*.py** files - Test modules for each algorithm

### Solver Implementations

- **CBCPLEX.py** - CPLEX solver (default)
- **CBGLPK.py** - GLPK solver (open source alternative)
- **CBGLPKOLD.py** - Legacy GLPK wrapper

### Multi-core support

- **CBMultiCore.py** - Parallel FBA/FVA using multiprocessing
- **CBMultiModel.py** - Batch model operations
- **CBMultiEnv.py** - Environment management for parallel runs

## Key Data Structures

### Metabolic Model (`Fbase`)
```python
class Fbase:
    - metabolites: list of metabolite objects
    - reactions: list of reaction objects
    - compartments: dictionary of compartments
    - gene_products: gene to metabolite associations
    - sbml_file: SBML file source
    - solution: LP solver results

    Methods:
    - doFBA() - Flux Balance Analysis
    - doFVA() - Flux Variability Analysis
    - getFVAdata() - Get flux bounds for a reaction
    - clone() - Clone model
    - getMetabolites() / getReactions()
```

### Reaction object
```python
- getId() - Reaction ID
- getLowerBound() / getUpperBound() - Flux bounds
- getFVAdata() - Returns [min_flux, max_flux, flux, span]
- toggleReversibility() - Enable/disable reversibility
```

### Metabolite object
```python
- getId() - Metabolite ID
- getCompartment() - Compartment name
- getGeneIds() - Associated gene products
```

## SBML Support

CBMPy supports:
- SBML Level 2 (FBA) - `readSBML2FBA()`
- SBML Level 3 FBC - `readSBML3FBC()`
- FBC V3 (quadratic objectives) - via CBSolverX

Write functions:
- `writeSBML3FBC()` - Standard FBC
- `writeSBML3FBCV3()` - FBC V3 with FBC annotations
- `writeSBML3FBCV2()` - FBC V2 (legacy)

## Sphinx Documentation

To build Sphinx documentation:

```bash
# Install sphinx dependencies
pip install -e ".[dev]"

# Generate API documentation from docstrings
sphinx-apidoc -o docs/source cbmpy/

# Build the documentation
cd docs
make html

# Open the generated docs in browser
xdg-open build/html/index.html
```

## License

GNU General Public License v3 or later (GPLv3+)

## Contact

- Developer: Brett G. Olivier (b.g.olivier@vu.nl)
- GitHub: https://github.com/SystemsBioinformatics/cbmpy
- Issue tracker: https://github.com/SystemsBioinformatics/cbmpy/issues
