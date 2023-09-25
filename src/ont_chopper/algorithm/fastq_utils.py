#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""FASTQ utils functions."""
from dataclasses import dataclass
from Bio import SeqIO
from array import array


@dataclass
class FastqRecord:
    """Class to store the information of FastqRecord."""

    uid: str
    name: str
    seq: str
    quality_str: str
    quality: array

    def __repr__(self) -> str:
        """Get a string representation of an FastqRecord."""
        return f"FastqRecord({self.name}:{self.seq})"

    def __hash__(self) -> int:
        """Hash an FastqRecord."""
        return hash(self.uid) ^ hash(self.seq) ^ hash(self.quality_str)

    def revcomp_seq(self):
        """Reverse complement sequence record."""
        return FastqRecord(
            self.uid,
            self.name,
            reverse_complement(self.seq),
            self.quality_str[::-1],
            self.quality[::-1],
        )

    def __len__(self) -> int:
        """Obtain reads length."""
        return len(self.seq)


def read_fq(fastq_file):
    """Read fastq files.

    This is a generator function that yields FastqRecord objects.
    """

    for record in SeqIO.parse(fastq_file, "fastq"):
        _, seq, _, quality_str, _ = record.format("fastq").split("\n")
        phred = array("i", record.letter_annotations["phred_quality"])
        seq_id = record.id
        seq_name = record.description
        yield FastqRecord(
            uid=seq_id,
            name=seq_name,
            seq=seq,
            quality_str=quality_str,
            quality=phred,
        )


def write_fq(read, fh):
    """Write read to fastq file."""
    fh.write(f"@{read.name}\n{read.seq}\n+\n{read.quality_str}\n")


def reverse_complement(seq: str) -> str:
    """Obtain reverse complement sequence."""
    rctrans = str.maketrans("ACGT", "TGCA")
    return str.translate(seq, rctrans)[::-1]
