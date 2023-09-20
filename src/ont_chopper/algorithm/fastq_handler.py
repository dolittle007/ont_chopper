#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================
import sys
import re
import os
from Bio import SeqIO
from .score_handler import valley_finder
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from array import array
from loguru import logger


def adapter_finder(
    record,
    score_cutoff,
    minimum_adapter_length,
    minimum_seq_length,
    minperc,
    window,
    queue,
):
    score = record.letter_annotations["phred_quality"]
    seq_info, seq, _, quality, _ = record.format("fastq").split("\n")
    read_length = len(seq)
    if " " in seq_info:
        _seq_id, description = seq_info.split(" ", 1)
    else:
        _seq_id = seq_info

    seq_id = _seq_id.lstrip("@")
    logger.info(f"{seq_id=}")
    adapter_positions = valley_finder(score, score_cutoff, minperc, window)

    if adapter_positions:
        position_list = array("l", [0])
        for valley_start, valley_end in adapter_positions:
            if valley_end - valley_start > minimum_adapter_length:
                position_list.append(valley_start)
                position_list.append(valley_end)
                logger.info(
                    f"adapter_start: {valley_start}, adapter_end: {valley_end}"
                )
        position_list.append(read_length - 1)

        if len(position_list) == 2:
            output_record = f"{seq_info}\n{seq}\n+\n{quality}\n"
        else:
            output_record = ''
            for segment_start, segment_end in zip(position_list[0::2], position_list[1::2]):
                if segment_end - segment_start > minimum_seq_length:
                    output_record += (
                        f"@{segment_start}:{segment_end}|{seq_info}\n"
                        f"{seq[segment_start:segment_end]}\n+\n"
                        f"{quality[segment_start:segment_end]}\n"
                    )
        del position_list
    else:
        output_record = f"{seq_info}\n{seq}\n+\n{quality}\n"

    if output_record:
        queue.put(output_record)
    return output_record


def listener(queue, out_filename):
    """listens for messages on the q, writes to file."""

    with open(out_filename, "w") as f:
        while True:
            contents = queue.get()
            if contents is None:
                break
            f.write(contents)
            f.flush()


def fastq_io(
    in_fastq: str,
    out_fastq: str,
    threads_num: int,
    minimum_seq_length: int = 150,
    minimum_adapter_length: int = 30,
    score_cutoff: int = 20,
    minperc: int = 10,
    window: int = 50,
) -> str:
    """ """
    cpu_number = multiprocessing.cpu_count()
    reads_number = 0
    found_adapter_number = 0
    effective_threads_num = min(threads_num, cpu_number)
    logger.info(f"{cpu_number=}, {effective_threads_num=}")

    if os.path.exists(out_fastq):
        os.remove(out_fastq)

    manager = multiprocessing.Manager()
    queue = manager.Queue()
    pool = multiprocessing.Pool(effective_threads_num)

    watcher = pool.apply_async(listener, (queue, out_fastq))

    jobs = []
    for record in SeqIO.parse(in_fastq, "fastq"):
        job = pool.apply_async(
            adapter_finder,
            (
                record,
                score_cutoff,
                minimum_adapter_length,
                minimum_seq_length,
                minperc,
                window,
                queue,
            ),
        )
        jobs.append(job)
        reads_number += 1
    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()
    # now we are done, kill the listener
    queue.put(None)
    pool.close()
    pool.join()

    logger.info(f"total reads number: {reads_number}")
    return out_fastq
