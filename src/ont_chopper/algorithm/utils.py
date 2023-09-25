#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""utils for ont-chopper."""
from itertools import islice, chain


def count_fastq_records(fname, size=128000000, opener=open):
    """Obtain number of records in FASTQ file."""
    fh = opener(fname, "r")
    count = 0
    while True:
        b = fh.read(size)
        if not b:
            break
        count += b.count("\n+\n")
    fh.close()
    return count


def batch_process(iterable, batch_size):
    """Process data in batch size."""
    sourceiter = iter(iterable)
    while True:
        batchiter = islice(sourceiter, batch_size)
        try:
            yield list(chain([next(batchiter)], batchiter))
        except StopIteration:
            return
