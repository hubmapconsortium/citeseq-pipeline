#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from operator import index
from os import environ, fspath
from pathlib import Path
from subprocess import check_call
from typing import Iterable, Optional, Sequence, Tuple

from fastq_utils import find_grouped_fastq_files

SALMON_COMMAND = [
    "salmon",
    "alevin",
    "-l",
    "ISR",
    "--index",
    "{index_dir}",
    "-o",
    "{out_dir}",
    "--citeseq",
    "--featureStart",
    "0",
    "--featureLength",
    "15",
    "-p",
    "{threads}",
    "--dumpMtx",
]


def main(
    fastq_dir: Path,
    threads: Optional[int],
    index_dir: Path,
):
    threads = threads or 1
    path = None
    print("Target index directory is", str(index_dir))
    if "adt_index" in str(index_dir):
        path = "alevin_adt"
        print("Starting ADT quantification...")
    elif "hto_index" in str(index_dir):
        path = "alevin_hto"
        print("Starting HTO quantification...")
    else:
        raise ValueError("Incorrect input index directory! Not for ADT/HTO experiments")

    command = [
        piece.format(
            salmon_option="citeseq",
            threads=threads,
            index_dir=index_dir,
            out_dir=path,
        )
        for piece in SALMON_COMMAND
    ]

    if "hto_index" in str(index_dir):
        command.extend(["--naiveEqclass"])
        print("Extended --naiveEqclass flag for HTO experiment.")

    fastq_pairs: Iterable[Sequence[Path]]

    fastq_pairs = list(find_grouped_fastq_files(fastq_dir, 2))

    if not fastq_pairs:
        raise ValueError("No FASTQ files found")

    for r1_fastq_file, r2_fastq_file in fastq_pairs:
        fastq_extension = [
            "-1",
            r1_fastq_file,
            "-2",
            r2_fastq_file,
        ]
        command.extend(fastq_extension)

    print("Running:", " ".join(str(x) for x in command))

    env = environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/local/lib"
    check_call(command)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--fastq_dir", type=Path)
    p.add_argument("-p", "--threads", type=int)
    p.add_argument("--index_dir", type=Path)
    args = p.parse_args()

    main(
        args.fastq_dir,
        args.threads,
        args.index_dir,
    )
