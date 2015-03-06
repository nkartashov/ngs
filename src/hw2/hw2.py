__author__ = 'nikita_kartashov'

from math import ceil
from sys import argv
import pysam
from os import path, makedirs

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
        print("[{0}, {1}]".format(left, right))
        for reference in bam.references:
            coverage[i] += bam.count(reference=reference, start=left, end=right)
        coverage[i] = coverage[i] * 1.0 / AVERAGE_STEP
    plt.plot(range(0, len(coverage)), coverage)
    plt.ylabel('Coverage')
    plt.xlabel('Base x1000')
    plt.axis([0, len(coverage), 0, max(coverage)])
    plt.savefig(result_file('base_coverage.png'))
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
            pass
