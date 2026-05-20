# Install CBMPy, with MiniConda and CPLEX

## Python environment
Install miniconda from: https://www.anaconda.com/docs/getting-started/miniconda/main

When using any CPLEX 22.1.2 or newer with DOCPLEX (recommended)
- `conda create -n cbmpy python=3.12 ipython pip wheel numpy scipy matplotlib`

Or alternatively when using CPLEX 22.1.1 (not recommended)
- `conda create -n cbmpy python=3.10 ipython pip wheel numpy scipy matplotlib`

Then activate your conda evironment with:
- `conda activate cbmpy`
- `pip install python-libsbml xlrd xlwt biopython xlsxwriter`

## Install Solver
GLPK will be installed but is limited in functionality, install CPLEX 22.1.2 (recommended) or use the CPLEX 22.1.11 with Python 3.10 workaround (not recommended)

### With CPLEX 22.1.2 and Python 3.12
Download and install: `cplex_studio2212.win_x86_64.exe` 
- https://academic.ibm.com/a2mt/downloads/data_science#/
- https://academic.ibm.com/a2mt/downloads/data_science?cardID=52

If necessary activate your CBMPy environment: `conda activate cbmpy` 
- `pip install cplex docplex`

Point docplex to your CPLEX studio 22.1.2 base directory
- `docplex config --upgrade c:\...\CPLEX_Studio2212\`

### Using old cplex 22.1.1 and Python 3.10 (not recommmended)
- `cd CPLEX_Studio2212\cplex\python\3.10\x64_win64`
- `pip install .`

## Install CBMPy
If necessary activate your CBMPy environment: `conda activate cbmpy` 

### With git installed
git clone https://github.com/SystemsBioinformatics/cbmpy.git
cd cbmpy
git checkout dev
git status
pip install .

### With source
Change directory to the CBMPy dev source and install:
```bash
cd cbmpy-dev-src
pip install .
```
