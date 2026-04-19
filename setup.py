"""
CBMPy: setup.py (legacy)

This file is kept for backward compatibility with legacy setup.py tools.
All build configuration is managed by pyproject.toml.

To build, use:
    python -m build
    pip install -e .
"""

from setuptools import setup

# Legacy: version is now read from pyproject.toml [project]
# This setup.py is intentionally minimal.

setup()
