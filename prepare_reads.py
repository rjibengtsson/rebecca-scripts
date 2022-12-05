'''
Description: This script is taken from parts of the customer pipeline to prepare fastq reads for downstream analysis.
It will perform adapters and quality trimming, deduplication and down sampling if coverage is above a threshold.

Usage: python prepare_reads.py -s1 ERR434473_1.fastq.gz -s2 ERR434473_2.fastq.gz -adapter bbmap-adapters.fa -codex ../../e_coli_v1.0.0.0.fasta -out ./processed/
'''

##################
# Import modules #
##################

import os
import gzip
import argparse
from pathlib import Path
import subprocess
import itertools
import typing as t

parser = argparse.ArgumentParser(
    prog = 'prepare_reads.py',
    usage='%(prog)s -s1 read_1.fastq.gz -s2 read_2.fastq.gz -adapter bbmap-adapters.fa -codex e_coli_v1.0.0.0.fasta -out OUT_DIR',
    description='Process fastq reads according to the customer pipeline')
parser.add_argument('-s1','--s1', help="read 1", required=True)
parser.add_argument('-s2','--s2', help="read 2", required=True)
parser.add_argument('-adapter','--adapter', help="adapter sequences", required=True)
parser.add_argument('-codex','--codex', help="codex file", required=True)
parser.add_argument('-out','--out', help="output directory", required=True)
args = parser.parse_args()

read_1 = args.s1
read_2 = args.s2
adapter_file = args.adapter
codex_file = args.codex
out_dir_path = args.out


# create directorys
out_dir = Path(out_dir_path)
os.mkdir(out_dir)
log_dir = Path(os.path.join(out_dir, "logs"))
os.mkdir(log_dir)


MINIMUM_COVERAGE = 5
MAXIMUM_COVERAGE = 200
TARGET_COVERAGE = 100
MIN_READ_COUNT = 100_000
BIG_GENOME = 7_000_000

sample_name = str(read_1).split("_")[0]

def extract_trimmed_stats(bbduk_log_path) -> dict[str, int]:
    stats = {}
    with open(bbduk_log_path, "r", encoding="utf8") as file:
        for line in file:
            if line.startswith("Result"):
                stats["bases_out"] = int(line.rstrip().split("\t")[2].split(" ")[0])
                stats["reads_out"] = int(line.rstrip().split("\t")[1].split(" ")[0])
            elif line.startswith("Input:"):
                stats["reads_in"] = int(line.rstrip().split("\t")[1].split(" ")[0])
                stats["bases_in"] = int(line.rstrip().split("\t")[3].split(" ")[0])

    missing_stats = {"bases_in", "bases_out", "reads_in", "reads_out"} - stats.keys()
    if missing_stats:
        raise Exception(f"Failed to extract bbduk stats: {missing_stats}")
    return stats


# Perform adapter trimming and quality filtering
def trim_readsets(read_1, read_2, adapter_file):
    trimmed_1 = out_dir.joinpath(f"{sample_name}_trimmed_1.fastq.gz")
    trimmed_2 = out_dir.joinpath(f"{sample_name}_trimmed_2.fastq.gz")
    cmd = [
        "bbduk.sh",
        "-Xmx10g",
        f"in1={read_1}",
        f"in2={read_2}",
        f"out1={trimmed_1}",
        f"out2={trimmed_2}",
        f"ref={adapter_file}",
        "ktrim=r",
        "k=23",
        "mink=11",
        "hdist=1",
        "tpe",
        "qtrim=rl",
        "trimq=20",
        "ftm=5",
        "minlen=72",
        "ordered",
        "overwrite=t",
    ]
    bbduk_log_path = log_dir.joinpath(f"{sample_name}_bbduk.log")
    with open(bbduk_log_path, "w", encoding="utf8") as log_file:
        subprocess.run(cmd, stderr=log_file, check=True, text=True)
    trimming_stats = extract_trimmed_stats(bbduk_log_path)
    return trimmed_1, trimmed_2, trimming_stats



CLUMPIFY_STATS = {
    "Bases Out": "bases_out",
    "Reads Out": "reads_out",
    "Reads In": "reads_in",
}


def extract_dedupe_stats(dedupe_log_path) -> dict[str, int]:
    stats = {}
    with open(dedupe_log_path, "r", encoding="utf8") as file:
        for line in file:
            for output_key, stat_key in CLUMPIFY_STATS.items():
                if line.startswith(f"{output_key}:"):
                    stats[stat_key] = int(line.rstrip().split(":")[1].lstrip(" "))
    missing_stats = set(CLUMPIFY_STATS.values()) - stats.keys()
    if missing_stats:
        raise Exception(f"Couldn't extract stats from deduplication: {missing_stats}")
    return stats


# Perform deduplication
def dedupe(trimmed_1, trimmed_2):
    dedupe_1 = out_dir.joinpath(f"{sample_name}_dedupe_1.fastq.gz")
    dedupe_2 = out_dir.joinpath(f"{sample_name}_dedupe_2.fastq.gz")
    cmd = [
        "clumpify.sh",
        "-Xmx4g",
        "t=2",
        f"in1={trimmed_1}",
        f"in2={trimmed_2}",
        f"out1={dedupe_1}",
        f"out2={dedupe_2}",
        "dedupe",
        "reorder",
    ]

    # print dedupe log
    dedupe_log_path = log_dir.joinpath(f"{sample_name}_dedupe.log")
    with open(dedupe_log_path, "w", encoding="utf8") as log_file:
        subprocess.run(cmd, stderr=log_file, check=True, text=True)
    stats = extract_dedupe_stats(dedupe_log_path)

    return dedupe_1, dedupe_2, stats["bases_out"]



# Perform downsampling
def downsample(dedupe_1, dedupe_2, size_or_fraction):
    downsample_1 = out_dir.joinpath(f"{sample_name}_downsample_1.fastq.gz")
    downsample_2 = out_dir.joinpath(f"{sample_name}_downsample_2.fastq.gz")
    cmd1 = [
        f"seqtk sample -s 1 {dedupe_1} {size_or_fraction} | gzip - > {downsample_1}"
    ]
    subprocess.run(cmd1, shell=True)
    cmd2 = [
        f"seqtk sample -s 1 {dedupe_2} {size_or_fraction} | gzip - > {downsample_2}"
    ]
    subprocess.run(cmd2, shell=True)
    return downsample_1, downsample_2


# Determine if reads need downsampling or not
def maybe_downsample(trimmed_1, trimmed_2):
    deduped_1, deduped_2, number_dedupe_bases = dedupe(trimmed_1, trimmed_2)
    coverage = number_dedupe_bases / 12542403
    multiplier = TARGET_COVERAGE / coverage
    if coverage > MAXIMUM_COVERAGE:
        return downsample(dedupe_1, deduped_2, multiplier)


# function to execute pipeline
def run(read_1, read_2, adapter_file) -> None:
    trimmed_1, trimmed_2, trimming_stats = trim_readsets(read_1, read_2, adapter_file)
    dedupe(trimmed_1, trimmed_2)
    maybe_downsample(trimmed_1, trimmed_2)

# execute pipeline
run(read_1, read_2, adapter_file)