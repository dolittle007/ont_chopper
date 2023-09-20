#!/usr/bin/env python
# -*- coding: utf-8 -*-
from findpeaks import findpeaks
from typing import List
from typing import Tuple
from typing import Union


def obatin_adjacent_peak_index(tgt_index, peak_candidates_indices, direction="greater"):
    sorted_peak_candidates_indeces = sorted(peak_candidates_indices)
    adjacent_peak_index = None
    if direction == "greater":
        for _idx in sorted_peak_candidates_indeces:
            if _idx > tgt_index:
                adjacent_peak_index = _idx
                break
    elif direction == "less":
        for _idx in sorted_peak_candidates_indeces[::-1]:
            if _idx < tgt_index:
                adjacent_peak_index = _idx
                break
    return adjacent_peak_index


def valley_extension(
    tgt_index, phred_quality_score, peak_candidates_indices, score_cutoff
):
    read_length = len(phred_quality_score)
    upper_bound = obatin_adjacent_peak_index(
        tgt_index, peak_candidates_indices, direction="greater"
    )
    lower_bound = obatin_adjacent_peak_index(
        tgt_index, peak_candidates_indices, direction="less"
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
