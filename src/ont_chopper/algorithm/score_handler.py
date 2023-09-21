#!/usr/bin/env python
# -*- coding: utf-8 -*-
from findpeaks import findpeaks
from typing import List
from typing import Tuple
from typing import Union


def find_numbers_less_and_greater(numbers, target) -> Tuple:
    """the first number that is less than and the first number that
        is greater than a target number in a list of numbers
    :param numbers: a list of numbers
    :param target: the target number
    :type numbers: list
    :type target: int
    :rtype: tuple
    """

    less_than_target = None
    greater_than_target = None

    for num in numbers:
        if num < target:
            if less_than_target is None or num > less_than_target:
                less_than_target = num
        elif num > target:
            if greater_than_target is None or num < greater_than_target:
                greater_than_target = num

    return less_than_target, greater_than_target


def valley_extension(
    tgt_index, phred_quality_score, peak_candidates_indices, score_cutoff
):
    """ Obtain valley zone that satisfy score cutoff.
    """
    read_length = len(phred_quality_score)

    lower_bound, upper_bound = find_numbers_less_and_greater(
        peak_candidates_indices, tgt_index
    )

    if upper_bound is None:
        upper_bound = read_length - 1

    if lower_bound is None:
        lower_bound = 0

    # extend in the right direction
    tgt_end = tgt_index
    for _idx in range(tgt_index, upper_bound + 1):
        if phred_quality_score[_idx] <= score_cutoff:
            tgt_end += 1
        else:
            break

    # extend in the left direction
    tgt_start = tgt_index
    for _idx in range(tgt_index, lower_bound - 1, -1):
        if phred_quality_score[_idx] <= score_cutoff:
            tgt_start -= 1
        else:
            break

    return (tgt_start, tgt_end)


def valley_finder(
    phred_quality_score: List[int],
    score_cutoff: int = 20,
    minperc: int = 10,
    window: int = 50,
) -> List[Tuple[int, int]]:
    """Identify valley in phred_quality_score."""

    adapter_positions = []
    if len(phred_quality_score) < window * 2:
        return adapter_positions

    fp = findpeaks(
        method="caerus", params_caerus={"minperc": minperc, "window": window}, verbose=0
    )
    result = fp.fit(phred_quality_score)
    data = result["df"]

    del fp
    del result

    valley_positions = list(
        data[(data["valley"] == True) & (data["y"] <= score_cutoff)]["x"]
    )
    peak_positions = list(
        data[(data["peak"] == True) & (data["y"] > score_cutoff)]["x"]
    )

    for tgt_index in valley_positions:
        position_tuple = valley_extension(
            tgt_index, phred_quality_score, peak_positions, score_cutoff
        )
        if not position_tuple in adapter_positions:
            adapter_positions.append(position_tuple)
    return adapter_positions
