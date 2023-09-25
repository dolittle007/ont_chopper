#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fastq IO module."""
from concurrent.futures import ProcessPoolExecutor
from array import array
from loguru import logger
import tqdm
from typing import List
from typing import Tuple
from .utils import count_fastq_records, batch_process
from .fastq_utils import read_fq, write_fq, FastqRecord
from .adapter_finder import adapter_search


def segments_to_reads(read, segments, minimum_seq_length):
    """Convert segments to output reads with annotation."""
    for start, end in segments:
        if end - start > minimum_seq_length:
            yield FastqRecord(
                uid=f"{start}:{end}|{read.uid}",
                name=read.name,
                seq=read.seq[start:end],
                quality_str=read.quality_str[start:end],
                quality=read.phred[start:end],
            )


def _find_adapters_single(params) -> List[Tuple[int, int]]:
    """Find non-adapter regions in a single reads."""
    read = params[0]
    (
        phred_threshold,
        minimum_adapter_len,
        maximum_adapter_len,
        required_polya_len,
    ) = params[1]

    adapters = adapter_search(
        read.seq,
        read.quality,
        required_polya_len,
        phred_threshold,
        minimum_adapter_len,
        maximum_adapter_len,
    )

    read_length = len(read)
    segments = []

    if len(adapters) > 0:
        _positions = [0]
        for i, j in adapters:
            if j != read_length:
                _positions.extend([i, j])
            else:
                _positions.append(i)
        segments = [(x, y) for x, y in zip(_positions[::2], _positions[1::2])]
    return segments


def find_adapters(
    reads,
    pool,
    min_batch_size,
    phred_threshold,
    minimum_adapter_len,
    maximum_adapter_len,
    required_polya_len,
):
    for batch in batch_process(reads, min_batch_size):
        for res in pool.map(
            _find_adapters_single,
            zip(
                batch,
                [
                    (
                        phred_threshold,
                        minimum_adapter_len,
                        maximum_adapter_len,
                        required_polya_len,
                    )
                ]
                * len(batch),
            ),
        ):
            try:
                yield res
            except StopIteration:
                return


def chopper(
    reads,
    pool,
    min_batch_size,
    phred_threshold,
    minimum_adapter_len,
    maximum_adapter_len,
    required_polya_len,
):
    batch_segments = find_adapters(
        reads,
        pool,
        min_batch_size,
        phred_threshold,
        minimum_adapter_len,
        maximum_adapter_len,
        required_polya_len,
    )
    for i, segments in enumerate(batch_segments):
        yield reads[i], segments


def fastq_io(
    in_fastq: str,
    rescued_fastq: str,
    unclassified_fastq: str,
    threads_num: int,
    batch_size: int,
    logger,
    minimum_seq_length: int = 150,
    phred_threshold: int = 20,
    minimum_adapter_len: int = 30,
    maximum_adapter_len: int = 100,
    required_polya_len: int = 2,
):
    """Read FASTQ and Write FASTQ files."""
    records_num = count_fastq_records(in_fastq)

    if batch_size > records_num:
        batch_size = records_num
    if batch_size == 0:
        batch_size = 1

    logger.info(f"Processing the whole dataset using a batch size of {batch_size}\n")

    progress_bar = tqdm.tqdm(total=records_num)
    min_batch_size = max(batch_size / threads_num, 1)

    unclassified_fh = open(unclassified_fastq, "w")
    rescued_fh = open(rescued_fastq, "w")
    with ProcessPoolExecutor(max_workers=threads_num) as executor:
        for batch in batch_process(read_fq(in_fastq), batch_size):
            for read, segments in chopper(batch, executor, min_batch_size):
                # reads without identified adapters
                if unclassified_fastq is not None and len(segments) == 0:
                    write_fq(read, unclassified_fh)
                # reads with identified adapters
                if rescued_fastq is not None and len(segments) > 0:
                    for trim_read in segments_to_reads(read, segments, minimum_seq_length):
                        write_fq(trim_read, rescued_fh)
                progress_bar.update(1)
    progress_bar.close()
    logger.info(f"Finished processing file: {in_fastq}\n")

    for fh in (
        unclassified_fh,
        rescued_fh,
    ):
        if fh is None:
            continue
        fh.flush()
        fh.close()
