# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# CBMPy - Constraint Based Modelling in Python

CBMPy is a Python platform for constraint-based modelling and analysis of biological systems. It implements Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), and other constraint-based methods using linear and mixed-integer programming.

## Installation

**Prefer conda environments.**

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

# Run a specific test file
pytest tests/test_specific.py

# Run a specific test function
pytest tests/test_specific.py::test_something

# Run tests with coverage
pytest --cov=cbmpy

# Exclude slow tests
pytest -m "not slow"

# Run only integration tests
pytest -m integration

# Test across Python versions via tox
tox
```

### Building and packaging

```bash
# Build source distribution and wheel
python setup.py sdist bdist_wheel

# Clean and rebuild
rm -rf build dist *.egg-info
python setup.py sdist bdist_wheel
```

### Code style and linting

The project uses pylint with custom disable comments for:
- C0103 (invalid variable name)
- C0301 (line too long)
- E1101 (module has no member)

Configured in `pyproject.toml` under `[tool.pylint.messages_control]`.

Docstrings use **numpydoc** style, configured in `pyproject.toml` under `[tool.docformatter]`.
Format all docstrings:

```bash
docformatter --in-place --recursive --style=numpy cbmpy/
```

Check formatting:

```bash
docformatter --diff --recursive cbmpy/
```

## Project Architecture

### Main Package: `cbmpy/`

**Core modules:**

- **CBModel.py** - Core model class. Contains `Fbase` base class for all model objects. Manages metabolites, reactions, and compartments. Has SBML read/write functionality.
- **CBSolver.py** - Solver orchestrator. Loads GLPK and/or CPLEX at import time (controlled by `CBMPY_USE_GLPK`/`CBMPY_USE_CPLX` env vars). Main entry points: `analyzeModel` (FBA), `FluxVariabilityAnalysis` (FVA).
- **CBConfig.py** - Single configuration dict (`__CBCONFIG__`) with version, solver prefs, and global constants.
- **CBRead.py / CBWrite.py** - SBML file I/O. `readSBML3FBC`, `readSBML2FBA`, `readCOBRASBML` for reading. `writeSBML3FBC`, `writeSBML3FBCV2`, `writeSBML3FBCV3`, `writeCOBRASBML` for writing.
- **CBCommon.py** - Common utilities: `StructMatrixLP`, MIRIAM annotation parsing, ID checking/fixing.
- **CBSolverX.py** - Extended solver with QP (quadratic objective) support and additional constraint types.
- **CBTools.py** - General utility functions (serialization, deSerialization).
- **CBXML.py** - SBML handling layer over python-libsbml.
- **CBNetDB.py** - Network analysis utilities.
- **CBPlot.py** - Visualization.
- **CBMultiCore.py / CBMultiModel.py / CBMultiEnv.py** - Multi-processing support for parallel FBA/FVA and batch model operations.
- **CBGUI.py / CBWx.py / CBQt4.py** - GUI frontends (wxWidgets and Qt4).

### Flux Modules: `cbmpy/fluxmodules/`

Advanced flux decomposition analyses:

- **fluxmodules.py** - Main entry point with `computeModules()` and `computeModulesMatroid()`. Requires prior FVA.
- **matroid.py** - Matroid theory implementation.
- **sparserationals.py** - Sparse rational number arithmetic.
- **decomposition.py / enumerate.py** - Flux module decomposition and enumeration algorithms.
- **random_color.py** - Graph coloring utilities.

### Solver Implementations

- **CBCPLEX.py** - IBM CPLEX solver (preferred, higher performance).
- **CBGLPK.py** - GLPK solver (open source default, loaded first by default).
- **CBGLPKOLD.py** - Legacy GLPK wrapper (kept for compatibility).

### Key Data Structures

- **`Fbase`** (in CBModel.py) - Base class for metabolic models. Tracks metabolites, reactions, compartments, gene products. Core methods: `doFBA()`, `doFVA()`, `getFVAdata()`, `clone()`.
- **Reaction** - Has `getId()`, `getLowerBound()`, `getUpperBound()`, `getFVAdata()`, `toggleReversibility()`.
- **Metabolite** - Has `getId()`, `getCompartment()`, `getGeneIds()`.

### Top-level API (`import cbmpy`)

The package `__init__.py` flattens the API:

```python
import cbmpy

# Reading models
cmod = cbmpy.readSBML3FBC('model.xml')

# Analyses
cbmpy.doFBA(cmod)      # Flux Balance Analysis
cbmpy.doFVA(cmod)      # Flux Variability Analysis

# Writing models
cbmpy.writeSBML3FBC(cmod, 'output.xml')
```

## Configuration

All global config lives in `__CBCONFIG__` (CBConfig.py). Key entries:
- `SOLVER_PREF` - Default solver preference (default: `'CPLEX'`)
- `SOLVER_ACTIVE` - Currently loaded solver(s)
- `DEBUG` - Debug flag
- `SYMPY_DENOM_LIMIT` - Symbolic computation denominator limit

Override solver via env vars: `CBMPY_USE_GLPK=1` or `CBMPY_USE_CPLX=1`.

## pyproject.toml Configuration

Build, test, and tooling config is all in `pyproject.toml`:
- `[tool.pytest.ini_options]` - Test paths (`tests/`), markers (`slow`, `integration`)
- `[tool.coverage.run]` / `[tool.coverage.report]` - Coverage settings
- `[tool.docformatter]` - numpydoc style, line width 79, black-compatible
- `[tool.pylint.messages_control]` - Disabled linter codes

## License

GNU General Public License v3 or later (GPLv3+)

## Contact

- Developer: Brett G. Olivier (@bgoli)
- GitHub: https://github.com/SystemsBioinformatics/cbmpy
- Issue tracker: https://github.com/SystemsBioinformatics/cbmpy/issues
