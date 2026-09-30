"""Argument parsing and dispatch for the crunchflow command-line interface."""

import argparse
from importlib.metadata import PackageNotFoundError, version

from crunchflow.util import SPATIAL_PROFILE_SUFFIXES, clear_output

try:
    __version__ = version("crunchflow")
except PackageNotFoundError:  # Running from a source tree rather than an install
    __version__ = "unknown"


def _add_clear_output_parser(subparsers):
    """Register the ``clear-output`` subcommand on `subparsers`."""
    parser = subparsers.add_parser(
        "clear-output",
        help="delete CrunchFlow output files from a run folder",
        description=(
            "Delete the output files that CrunchFlow writes into a run folder, so that a "
            "simulation can be re-run from a clean directory. Files are matched by name: a "
            "file is deleted just when its name is a known CrunchFlow output prefix followed "
            "by a time-step index and a recognized suffix. Subdirectories are left untouched, "
            "and time_series files are not removed, because CrunchFlow takes their names from "
            "the input file rather than generating them."
        ),
    )
    parser.add_argument(
        "-f",
        "--folder",
        default=".",
        help="folder to clear (default: the current directory)",
    )
    parser.add_argument(
        "-s",
        "--suffixes",
        nargs="+",
        metavar="SUFFIX",
        default=list(SPATIAL_PROFILE_SUFFIXES),
        help="output file suffixes to delete (default: %(default)s)",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="list the files that would be deleted, without deleting them",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print the list of deleted files",
    )
    parser.set_defaults(func=_run_clear_output)
    return parser


def _run_clear_output(args):
    """Run the ``clear-output`` subcommand and return its exit status."""
    deleted = clear_output(
        folder=args.folder,
        suffixes=args.suffixes,
        dry_run=args.dry_run,
        verbose=False,
    )

    if args.verbose:
        for path in deleted:
            print("    {}".format(path))
        verb = "Would delete" if args.dry_run else "Deleted"
        print("{} {} output file(s) in {}".format(verb, len(deleted), args.folder))

    return 0


def build_parser():
    """Build the top-level crunchflow argument parser.

    Returns
    -------
    argparse.ArgumentParser
        parser with every crunchflow subcommand registered on it
    """
    parser = argparse.ArgumentParser(
        prog="python -m crunchflow.cli",
        description="Command-line tools for the CrunchFlow reactive transport code.",
    )
    parser.add_argument("--version", action="version", version="crunchflow {}".format(__version__))

    subparsers = parser.add_subparsers(dest="command", metavar="command")
    _add_clear_output_parser(subparsers)

    return parser


def main(argv=None):
    """Run the crunchflow command-line interface.

    Parameters
    ----------
    argv : list of str, optional
        arguments to parse. The default is to read them from ``sys.argv``

    Returns
    -------
    int
        process exit status; 0 on success
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # No subcommand given, so there is nothing to dispatch to
    if getattr(args, "func", None) is None:
        parser.print_help()
        return 1

    return args.func(args)