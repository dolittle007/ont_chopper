# !/usr/bin/env python
"""CLi for ont_chopper.
"""
import argparse
import sys
from pathlib import Path
from typing import Union

from loguru import logger

from ..algorithm.fastq_io import fastq_io
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

    fastq_io(
        options.input_fastq,
        options.rescue_output,
        options.unclass_output,
        options.thread_num,
        options.batch_size,
        logger,
        options.minimum_seq_length,
        options.phred_threshold,
        options.minimum_adapter_length,
        options.maximum_adapter_length,
        options.required_polya_len,
    )
