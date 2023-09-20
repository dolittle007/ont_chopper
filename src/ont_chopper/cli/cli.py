# !/usr/bin/env python
"""CLi for ont_chopper.
"""
import argparse
import sys
import tempfile
import time
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Union

from loguru import logger

from ..algorithm.fastq_handler import fastq_io
from .arg import DefaultOptions


def cli(options: Union[argparse.Namespace, DefaultOptions]):
    """Cli function."""
    # add logger
    logger.remove()
    logger.add(
        sys.stdout,
        level=options.log.upper(),
        format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
        "<level>{level: <8}</level> | "
        "<cyan>{line}</cyan> - <level>{message}</level>",
        enqueue=True,
        colorize=True,
        backtrace=False,
        diagnose=True,
    )
    fastq_io(options.input, options.output, options.thread, options.minimum_seq_length, options.minimum_adapter_length, options.score, options.minperc, options.window)

