# !/usr/bin/env python
"""Parse command line arguments."""

import argparse
import textwrap
from dataclasses import dataclass
from typing import Any
from typing import Optional
from typing import Tuple

from ont_chopper import __version__


@dataclass
class DefaultOptions:
    """Cli default options."""

    input: str
    score: int = 20
    thread: int = 1
    minimum_seq_length: int = 150
    minimum_adapter_length: int = 30
    minperc: int = 10
    window: int = 50
    log: str = "info"
    output: str = "result.fastq"


class RichArgParser(argparse.ArgumentParser):
    """RichArgParser."""

    def __init__(self, *args: Any, **kwargs: Any):
        """RichArgParser."""
        from rich.console import Console

        self.console = Console()
        super().__init__(*args, **kwargs)

    @staticmethod
    def _color_message(message: str, color: str = "green") -> str:
        """Color message."""
        import re

        pattern = re.compile(r"(?P<arg>-{1,2}[-|\w]+)")
        return pattern.sub(lambda m: f"[bold {color}]{m.group('arg')}[/]", message)

    def _print_message(self, message: Optional[str], file: Any = None) -> None:
        if message:
            self.console.print(self._color_message(message))


class RichHelpFormatter(argparse.HelpFormatter):
    """RichHelpFormatter."""

    def __init__(self, *args: Any, **kwargs: Any):
        """RichHelpFormatter."""
        super().__init__(*args, max_help_position=42, **kwargs)  # type: ignore


def parse_args() -> argparse.ArgumentParser:
    """Parse command line arguments."""
    parser = RichArgParser(
        description="[red]ont_chopper[/] :rocket: "
        "Nanopore directRNA adapter remover",
        epilog=textwrap.dedent(
            """Authors: TingYou Wang,
            Northwestern University, 2022"""
        ),
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input FASTQ file",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--score",
        action="store",
        dest="score",
        help="maximum phred score of adapter sequence (default: %(default)s)",
        type=int,
        default=DefaultOptions.score,
    )
    parser.add_argument(
        "-t",
        "--thread",
        action="store",
        dest="thread",
        help="thread number (default: %(default)s)",
        type=int,
        default=DefaultOptions.thread,
    )
    parser.add_argument(
        "--min-seq-len",
        type=int,
        dest="minimum_seq_length",
        default=DefaultOptions.minimum_seq_length,
        help=f"minimum length of FASTQ sequence after chopping (default: %(default)s).",
    )
    parser.add_argument(
        "--min-adapter-len",
        type=int,
        dest="minimum_adapter_length",
        default=DefaultOptions.minimum_adapter_length,
        help=f"minimum length of adapter sequence (default: %(default)s).",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="Output adapter removed FASTQ file. (default: %(default)s)",
        default=DefaultOptions.output,
    )
    parser.add_argument(
        "--minperc",
        action="store",
        dest="minperc",
        type=int,
        help="peaks vallyes with at least percentage difference (default: %(default)s)",
        default=DefaultOptions.minperc,
    )
    parser.add_argument(
        "--window",
        action="store",
        dest="window",
        type=int,
        help="window size (default: %(default)s)",
        default=DefaultOptions.window,
    )
    parser.add_argument(
        "--log-level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace"],  # "warning", "error", "critical"
        default=DefaultOptions.log,
        help="set log level (default: %(default)s)",
    )

    parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    return parser
