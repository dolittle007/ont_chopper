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

    input_fastq: str
    unclass_output: str
    rescue_output: str
    phred_threshold: int = 20
    thread_num: int = 1
    batch_size: int = 10000
    minimum_seq_length: int = 150
    minimum_adapter_length: int = 30
    maximum_adapter_length: int = 100
    required_polya_len: int = 2
    log: str = "info"


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
        "Tool to identify and rescue full-length Nanopore direct RNA reads",
        epilog=textwrap.dedent(
            """Authors: TingYou Wang,
            Northwestern University, 2023"""
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
        "-u",
        "--unclass-output",
        action="store",
        dest="unclass_output",
        help="Write unclassified reads to this file. (Cannot find adapters)",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--rescue-output",
        action="store",
        dest="rescue_output",
        help="Write rescued reads to this file",
        required=True,
    )
    parser.add_argument(
        "-q",
        "--quality",
        action="store",
        dest="phred_threshold",
        help="Phred score cutoff of adapter sequence (default: %(default)s)",
        type=int,
        default=DefaultOptions.phred_threshold,
    )
    parser.add_argument(
        "-t",
        "--thread",
        action="store",
        dest="thread_num",
        help="thread number (default: %(default)s)",
        type=int,
        default=DefaultOptions.thread_num,
    )

    parser.add_argument(
        "-b",
        "--batch",
        action="store",
        dest="batch_size",
        help="Maximum number of reads processed in each batch (default: %(default)s)",
        type=int,
        default=DefaultOptions.batch_size,
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
        "--max-adapter-len",
        type=int,
        dest="maximum_adapter_length",
        default=DefaultOptions.maximum_adapter_length,
        help=f"maximum length of adapter sequence (default: %(default)s).",
    )

    parser.add_argument(
        "--min-polya-length",
        action="store",
        dest="required_polya_len",
        type=int,
        help="Minimum poly(A) tail length (default: %(default)s)",
        default=DefaultOptions.required_polya_len,
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
