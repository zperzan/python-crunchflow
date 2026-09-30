"""
Command-line interface to the crunchflow package.

Bundles the package's command-line tools as subcommands of a single parser,
invoked as ``python -m crunchflow.cli``. Run ``python -m crunchflow.cli --help``
for the list of available subcommands.
"""

from crunchflow.cli.main import main

__all__ = ["main"]
