#!/usr/bin/env python3
"""
Script to reformat docstrings to numpydoc format compatible with Sphinx.

This script uses docformatter to format docstrings to a consistent numpydoc style.
"""

import os
import re
import sys
from pathlib import Path

try:
    import tomli
except ImportError:
    import tomli_w
    import tomllib as tomli

try:
    import docformatter
except ImportError:
    from docformatter import __version__ as doc_version
    print(f"docformatter version: {doc_version}")


def get_project_config():
    """Get configuration from pyproject.toml or use defaults."""
    config_path = Path(__file__).parent / 'pyproject.toml'

    if config_path.exists():
        with open(config_path, 'rb') as f:
            config = tomli.load(f)
        docformatter_config = config.get('tool', {}).get('docformatter', {})
    else:
        docformatter_config = {}

    return docformatter_config


def format_module_file(filepath, exclude_patterns=None):
    """
    Format docstrings in a single file using docformatter.

    Args:
        filepath: Path to the Python file to format
        exclude_patterns: List of patterns to exclude from formatting

    Returns:
        True if file was formatted, False otherwise
    """
    try:
        if exclude_patterns:
            for pattern in exclude_patterns:
                if pattern in filepath:
                    print(f"Skipping {filepath} (matches exclusion pattern: {pattern})")
                    return False

        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        config = get_project_config()

        # Use docformatter to format the file
        result = docformatter.format_code(
            content,
            line_length=79,
            preselect_styles=['google', 'numpy', 'sphinx'],
            **config
        )

        if result != content:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(result)
            print(f"Formatted: {filepath}")
            return True

        return False

    except Exception as e:
        print(f"Error formatting {filepath}: {e}")
        return False


def find_all_python_files(root_path):
    """Find all Python files in the given path."""
    python_files = []
    exclude_patterns = [
        '__pycache__',
        '.git',
        '.tox',
        'build',
        'dist',
        '.pytest_cache',
        '.eggs',
    ]

    for filepath in Path(root_path).rglob('*.py'):
        filename = str(filepath)
        exclude = any(pattern in filename for pattern in exclude_patterns)
        if not exclude:
            python_files.append(filepath)

    return sorted(python_files)


def main():
    """Main function to format all docstrings."""
    root_path = Path(__file__).parent
    print(f"Scanning for Python files in {root_path}")

    files = find_all_python_files(root_path)
    print(f"Found {len(files)} Python files")

    formatted_count = 0

    for filepath in files:
        success = format_module_file(filepath)
        if success:
            formatted_count += 1

    print(f"\nFormatted {formatted_count} files")

    # Check for errors
    if formatted_count == 0:
        print("No files were formatted (may already be formatted or no docstrings found)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
