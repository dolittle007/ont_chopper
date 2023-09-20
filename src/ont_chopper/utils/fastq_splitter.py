#!/home/twp7981/Apps/Python3/envs/bio/bin/python
# -*- coding: utf-8 -*-
# =============================
import sys
import argparse
import os
from Bio import SeqIO
import shutil


def fastq_splitter(
    in_fastq: str,
    out_fastq_prefix: str,
    record_num: int = 1000000,
    temp_dir: str = "temp",
) -> str:
    """ """
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    else:
        os.mkdir(temp_dir)

    part_file = None
    part_num = 0
    record_index = 0
    for record in SeqIO.parse(in_fastq, "fastq"):
        if record_index % record_num == 0:
            part_num += 1
            if part_file:
                part_file.close()
            part_file_name = os.path.join(
                temp_dir, f"{out_fastq_prefix}_{part_num}.fastq"
            )
            part_file = open(part_file_name, "w")
        SeqIO.write(record, part_file, "fastq")
        record_index += 1
    if part_file:
        part_file.close()

    return True

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="input FASTQ file",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--number",
        action="store",
        dest="number",
        type=int,
        help="(default: %(default)s)",
        default=1000000,
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        action="store",
        dest="temp_dir",
        help="temp directory (default: %(default)s)",
        default="temp",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output_prefix",
        help="(default: %(default)s)",
        default="splitted",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()

    in_fastq = args.input
    out_fastq_prefix = args.output_prefix
    temp_dir = args.temp_dir
    record_num = args.number

    fastq_splitter(in_fastq, out_fastq_prefix, record_num, temp_dir)

if __name__ == "__main__":
    main()
