#!/usr/bin/env python
# -*- coding: utf-8 -*-
from array import array
import re
from statistics import median, mean
from loguru import logger
from typing import List
from typing import Tuple


def terminal_adapter_search(
    read_sequence: str,
    phred_quality_score: List[int],
    required_polya_len: int = 2,
    phred_threshold: int = 20,
    maximum_adapter_len: int = 100,
) -> Tuple[int, int]:
    """Identify putative adapter at 3' end, the average phred threshold applies for the sequence after polyA tail

    :param read_sequence: Read sequence
    :param phred_quality_score: Phred quality score
    :param required_polya_len: minimum required polyA length
    :param phred_threshold: phred score threshold
    :param maximum_adapter_len: maximum adapter length
    :type read_sequence: str
    :type phred_quality_score: array
    :type required_polya_len: int
    :type phred_threshold: int
    :type maximum_adapter_len: int
    :rtype: tuple

    adapters: [start, end)
    """
    read_length = len(read_sequence)
    polya_matches = [
        (match.start(), match.end())
        for match in re.finditer(rf"A{{{required_polya_len},}}", read_sequence)
    ]

    adapter_position_dict = {}
    for start, end in polya_matches:
        if (
            read_length - start <= maximum_adapter_len
        ):
            # read ends with polyA
            if end == read_length:
                adapter_position_dict[start] = end - start
            elif mean(phred_quality_score[end:]) <= phred_threshold:
                adapter_position_dict[start] = end - start


    # sorted by length of polyA length
    sorted_adapter_position_list = list(
        sorted(adapter_position_dict.items(), key=lambda item: item[1], reverse=True)
    )

    if len(sorted_adapter_position_list) > 0:
        adapter_position = sorted_adapter_position_list[0][0], read_length
    else:
        adapter_position = None, None

    return adapter_position


def internal_adapter_search(
    phred_quality_score: List[int],
    phred_threshold: int = 20,
    minimum_adapter_len: int = 30,
    distance_threshold: int = 1,
) -> List[Tuple[int, int]]:
    """Identify internal adapters in phred_quality_score.

    :param phred_quality_score: phred quality scores
    :param phred_threshold: phred score threshold
    :type phred_quality_score: array
    :type phred_threshold: int
    :rtype list

    adapters: [start, end)
    """

    reads_length = len(phred_quality_score)

    ascending_points = []
    descending_points = []

    for i in range(
        0, reads_length - 1
    ):  # from the first point to the second to last point
        if (
            phred_quality_score[i] <= phred_threshold
            and phred_quality_score[i + 1] > phred_threshold
        ):
            ascending_points.append(i)
        if (
            phred_quality_score[i] > phred_threshold
            and phred_quality_score[i + 1] <= phred_threshold
        ):
            descending_points.append(i + 1)
    # 5'end /
    if ascending_points[0] <= descending_points[0]:
        # 3' end /
        if ascending_points[-1] >= descending_points[-1]:
            descending_points.insert(0, 0)
        # 3' end \
        elif ascending_points[-1] < descending_points[-1]:
            descending_points.insert(0, 0)
            ascending_points.append(reads_length - 1)
    # 5' end \
    else:
        # 3' end /
        if ascending_points[-1] >= descending_points[-1]:
            pass
        # 3' end \
        elif ascending_points[-1] < descending_points[-1]:
            ascending_points.append(reads_length - 1)

    inital_adapter_positions = [
        (_start, _end) for _start, _end in zip(descending_points, ascending_points)
    ]

    merged_adapter_positions = merge_adjacent_regions(
        inital_adapter_positions, distance_threshold
    )

    adapter_positions = []
    for _start, _end in merged_adapter_positions:
        logger.debug(f"{_start=}, {_end=}")
        if _end - _start + 1 >= minimum_adapter_len:
            adapter_positions.append((_start, _end + 1))

    return adapter_positions


def merge_adjacent_regions(
    regions: List[Tuple[int, int]], distance_threshold: int = 1
) -> List[Tuple[int, int]]:
    """Merge genomic regions that are close to each other within a specified distance threshold.

    :param regions: putative adapter regions
    :param distance_threshold: distance threshold to define adjacent region
    :type regions: list
    :type distance_threshold: int
    :rtype list
    """

    # Initialize the merged regions list with the first region
    merged_regions = [regions[0]]

    # Iterate through the sorted regions
    for region in regions[1:]:
        prev_region = merged_regions[-1]  # Get the last merged region

        # Check if the current region is within the distance threshold of the previous one
        if region[0] - prev_region[1] <= distance_threshold + 1:
            # Merge the current region with the previous one
            merged_regions[-1] = (prev_region[0], max(prev_region[1], region[1]))
        else:
            # Add the current region as a new merged region
            merged_regions.append(region)

    return merged_regions


def adapter_search(
    read_sequence: str,
    phred_quality_score: List[int],
    required_polya_len: int = 2,
    phred_threshold: int = 20,
    minimum_adapter_len: int = 30,
    maximum_adapter_len: int = 100,
):
    adapters = []
    _terminal_adapter = terminal_adapter_search(
        read_sequence,
        phred_quality_score,
        required_polya_len,
        phred_threshold,
        maximum_adapter_len,
    )
    _internal_adapters = internal_adapter_search(
        phred_quality_score,
        phred_threshold,
        minimum_adapter_len,
    )
    if _terminal_adapter != (None, None):
        internal_adapters = [
            (i, j)
            for i, j in _internal_adapters
            if i <= _terminal_adapter[0] and j <= _terminal_adapter[1]
        ]
        if len(internal_adapters) > 0:
            adapters.extend(internal_adapters)
        adapters.append(_terminal_adapter)
    else:
        if len(_internal_adapters) > 0:
            adapters.extend(_internal_adapters)

    return adapters
