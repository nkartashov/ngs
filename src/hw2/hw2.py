__author__ = 'nikita_kartashov'

from math import ceil
from sys import argv
import pysam
from os import path, makedirs

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math

from Bio import SeqIO


def relative_path(relative_path_part):
    script_dir = path.dirname(__file__)
    return path.join(script_dir, relative_path_part)


RESULT_PATH = relative_path('../../results/hw2')


def result_file(filename):
    return path.join(RESULT_PATH, filename)


def make_result_dir():
    if not path.exists(RESULT_PATH):
        makedirs(RESULT_PATH)


def get_reference_length(sam):
    return sam.header["SQ"][0]["LN"]


AVERAGE_STEP = 1000


def get_genome_coverage(bam):
    length = get_reference_length(bam)
    regions = int(ceil(length * 1.0 / AVERAGE_STEP))
    coverage = [0 for _ in xrange(regions)]
    covered_regions = 0
    for i in xrange(len(coverage)):
        left = i * AVERAGE_STEP
        right = left + AVERAGE_STEP
        for reference in bam.references:
            coverage[i] += bam.count(reference=reference, start=left, end=right)
        coverage[i] = coverage[i] * 1.0 / AVERAGE_STEP
        if coverage[i] > 0:
            covered_regions += 1
    print("Covered genome rate: {0}".format(covered_regions * AVERAGE_STEP * 1.0 / length))
    plt.plot(range(0, len(coverage)), coverage)
    plt.ylabel('Coverage')
    plt.xlabel('Base x1000')
    plt.axis([0, len(coverage), 0, max(coverage)])
    plt.savefig(result_file('base_coverage.png'))
    plt.clf()


MAX_INSERTION_SIZE = 1000


def weighted_sd(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)  # Fast and numerically precise
    return math.sqrt(variance)


def get_insertion_size(bam):
    insertion_sizes = [0 for _ in xrange(MAX_INSERTION_SIZE)]

    for read in bam.fetch():
        if read.is_proper_pair and read.is_read1:
            insertion_sizes[abs(read.tlen)] += 1

    real_length = 0
    for i, size in enumerate(insertion_sizes):
        if size != 0:
            real_length = max(real_length, i)
    insertion_sizes = insertion_sizes[:real_length]
    sizes = [(i, number) for i, number in enumerate(insertion_sizes) if number != 0]
    sizes, weights = zip(*sizes)
    mean = np.average(sizes, weights=weights)
    sd = weighted_sd(sizes, weights=weights)
    quantiles = stats.norm.interval(0.95, loc=mean, scale=sd)
    print("IS mean: {0}\nIS sd: {1}".format(mean, sd))
    max_insertion_size = max(insertion_sizes)
    plt.plot(range(0, len(insertion_sizes)), insertion_sizes)
    plt.plot([mean, mean], [0, max_insertion_size], 'r-')
    for quantile in quantiles:
        plt.plot([quantile, quantile], [0, max_insertion_size], 'k-')
    plt.ylabel('Reads')
    plt.xlabel('Insertion size, base pairs')
    plt.axis([0, len(insertion_sizes), 0, max_insertion_size])
    plt.savefig(result_file('insertion_size.png'))
    plt.clf()


def get_mismatch_frequencies(bam):
    mismatches = {"A": {"C": 0, "T": 0, "G": 0},
                  "C": {"A": 0, "T": 0, "G": 0},
                  "G": {"A": 0, "C": 0, "T": 0},
                  "T": {"A": 0, "C": 0, "G": 0}}

    reference_genome = SeqIO.parse(open("../../resource/hw2/MG1655-K12.fasta"), 'fasta')
    genome = ""
    for fasta in reference_genome:
        genome = str(fasta.seq)

    for read in bam.fetch():
        reference_sequence = ""
        for start, end in read.get_blocks():
            reference_sequence += genome[start: end]

        read_sequence = read.query_alignment_sequence
        for i in xrange(len(reference_sequence)):
            if read_sequence[i] != reference_sequence[i] and \
               reference_sequence[i] != "N" and \
               read_sequence[i] != "N":
                mismatches[reference_sequence[i]][read_sequence[i]] += 1

    print("*\t{0}\t{1}\t{2}".format(mismatches["A"]["C"], mismatches["A"]["T"], mismatches["A"]["G"]))
    print("{0}\t*\t{1}\t{2}".format(mismatches["C"]["A"], mismatches["C"]["T"], mismatches["C"]["G"]))
    print("{0}\t{1}\t*\t{2}".format(mismatches["G"]["A"], mismatches["G"]["C"], mismatches["G"]["T"]))
    print("{0}\t{1}\t{2}\t*".format(mismatches["T"]["A"], mismatches["T"]["C"], mismatches["T"]["G"]))


if __name__ == '__main__':
    if len(argv) < 3:
        print("No paths to BAM files")
        exit(1)

    make_result_dir()

    bam_path = argv[1]
    mode = int(argv[2])
    tasks = {1: get_genome_coverage,
             2: get_insertion_size,
             3: get_mismatch_frequencies}
    if 1 > mode >= 3:
        print("Invalid task (valid 1-3)")
        exit(2)
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        tasks[mode](bamfile)