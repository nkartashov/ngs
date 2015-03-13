__author__ = 'nikita_kartashov'

from sys import argv
from os import makedirs, path
from math import ceil

import pysam
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from Bio import SeqIO

BWA_PREFIX = "bwa_"
BOWTIE_PREFIX = "bwt_"


def relative_path(relative_path_part):
    script_dir = path.dirname(__file__)
    return path.join(script_dir, relative_path_part)


RESULT_PATH = relative_path('../../results/hw3')


def result_file(filename):
    return path.join(RESULT_PATH, filename)


def make_result_dir():
    if not path.exists(RESULT_PATH):
        makedirs(RESULT_PATH)


def get_reference_length(sam):
    return sam.header["SQ"][0]["LN"]


AVERAGE_STEP = 1000


def get_genome_coverage(bam, prefix):
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
    mean = np.mean(coverage)
    print("Covered genome rate: {0}".format(covered_regions * AVERAGE_STEP * 1.0 / length))
    plt.plot(range(0, len(coverage)), coverage)
    plt.plot(range(0, len(coverage)), [mean for _ in xrange(0, len(coverage))], 'r-')
    plt.ylabel('Coverage')
    plt.xlabel('Base x1000')
    plt.axis([0, len(coverage), 0, max(coverage)])
    plt.savefig(result_file(prefix + 'base_coverage.png'))
    plt.clf()


def get_info_from_bam_reference(bam, reference_file, prefix=BWA_PREFIX):
    reference = SeqIO.parse(reference_file, 'fasta')
    get_genome_coverage(bam, prefix)

    pass


if __name__ == '__main__':
    if len(argv) < 4:
        print("No paths to files or mode")
        exit(1)

    make_result_dir()

    reference_path = argv[1]
    bamfile_path = argv[2]
    mode = int(argv[3])

    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile:
        prefix = ""
        if mode == 1:
            prefix = BWA_PREFIX
        elif mode == 2:
            prefix = BOWTIE_PREFIX
        else:
            print("Unrecognized mode")
            exit(2)
        with open(reference_path) as reference_file:
            get_info_from_bam_reference(bamfile, reference_file, prefix)
