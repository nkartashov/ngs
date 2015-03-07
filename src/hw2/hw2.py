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
    for i in xrange(len(coverage)):
        left = i * AVERAGE_STEP
        right = left + AVERAGE_STEP
        for reference in bam.references:
            coverage[i] += bam.count(reference=reference, start=left, end=right)
        coverage[i] = coverage[i] * 1.0 / AVERAGE_STEP
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


if __name__ == '__main__':
    if len(argv) < 3:
        print("No paths to BAM files")
        exit(1)

    make_result_dir()

    bam_path = argv[1]
    mode = argv[2] == "1"
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        if mode:
            get_genome_coverage(bamfile)
        else:
            get_insertion_size(bamfile)
