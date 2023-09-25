#!/usr/bin/env python3
"""Main function for ont chopper."""
import sys

from .cli.arg import parse_args
from .cli.cli import cli


def main():
    """Main function for dallas."""
    if sys.version_info < (3, 9):
        raise SystemExit(
            "Sorry, this code need Python 3.9 or higher. Please update. Aborting..."
        )

    parser = parse_args()
    options = parser.parse_args()

    cli(options)


if __name__ == "__main__":
    main()
