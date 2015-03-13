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


MAX_INDEL = 200

MATCH = 0
INSERTION = 1
DELETION = 2
SOFT_CLIPPED = 4
MISMATCH = 8


def get_indel_distribution(bam, prefix):
    indel_lengths = [0 for _ in xrange(MAX_INDEL)]
    max_indel = 0

    for read in bam.fetch():
        for t, length in read.cigartuples:
            if t == INSERTION or t == DELETION:
                indel_lengths[length] += 1
                max_indel = max(max_indel, length)

    indel_lengths = indel_lengths[1:max_indel]
    plt.plot(range(1, max_indel), indel_lengths)
    plt.ylabel('Indel length')
    plt.xlabel('Reads')
    plt.axis([1, max_indel, 0, max(indel_lengths)])
    plt.savefig(result_file(prefix + 'indel_distribution.png'))
    plt.clf()


def get_info_from_bam_reference(bam, reference_file, prefix=BWA_PREFIX):
    reference = (list(SeqIO.parse(reference_file, 'fasta'))[0]).seq
    get_genome_coverage(bam, prefix)
    # get homopolymer stuff
    get_indel_distribution(bam, prefix)
    get_quality_insertions(bam, prefix)
    get_quality_mismatches(bam, reference, prefix)
    get_mismatch_frequencies(bam, reference, prefix)


MAX_QUALITY = 50


def get_quality_insertions(bam, prefix):
    insertions_quality = [0 for _ in xrange(MAX_QUALITY)]
    max_quality = 0
    for read in bam.fetch():
        index = 0
        qualities = read.query_qualities
        for t, length in read.cigartuples:
            if t == INSERTION:
                for real_index in xrange(index, index + length):
                    insertions_quality[qualities[real_index]] += 1
                    max_quality = max(max_quality, qualities[real_index])
            if t == DELETION:
                continue
            index += length
    insertions_quality = insertions_quality[:max_quality]
    plt.plot(range(0, max_quality), insertions_quality)
    plt.ylabel('Insertions')
    plt.xlabel('Quality')
    plt.axis([0, max_quality, 0, max(insertions_quality)])
    plt.savefig(result_file(prefix + 'insertion_quality_distribution.png'))
    plt.clf()


def get_quality_mismatches(bam, reference, prefix):
    mismatch_quality = [0 for _ in xrange(MAX_QUALITY)]
    max_quality = 0
    for read in bam.fetch():
        reference_sequence = ""
        for start, end in read.get_blocks():
            reference_sequence += reference[start: end]
        read_sequence = read.query_alignment_sequence
        read_quality = read.query_alignment_qualities
        for i in xrange(len(reference_sequence)):
            if read_sequence[i] != reference_sequence[i] and \
                            reference_sequence[i] != "N" and \
                            read_sequence[i] != "N":
                mismatch_quality[read_quality[i]] += 1
                max_quality = max(max_quality, read_quality[i])
    mismatch_quality = mismatch_quality[:max_quality]
    plt.plot(range(0, max_quality), mismatch_quality)
    plt.ylabel('Mismatches')
    plt.xlabel('Quality')
    plt.axis([0, max_quality, 0, max(mismatch_quality)])
    plt.savefig(result_file(prefix + 'mismatch_quality_distribution.png'))
    plt.clf()


NUCLEOTIDES = {'A', 'C', 'T', 'G'}


def get_mismatch_frequencies(bam, reference_sequence, prefix):
    mismatches = {"A": {"C": 0, "T": 0, "G": 0, None: 0},
                  "C": {"A": 0, "T": 0, "G": 0, None: 0},
                  "G": {"A": 0, "C": 0, "T": 0, None: 0},
                  "T": {"A": 0, "C": 0, "G": 0, None: 0},
                  None: {"A": 0, "C": 0, "G": 0, "T": 0}}

    def filter_soft_clipped(cigar_tuples):
        return list(filter(lambda (t, l): t != SOFT_CLIPPED, cigar_tuples))

    for read in bam.fetch():
        reference_index = read.reference_start
        read_index = 0
        read_sequence = read.query_alignment_sequence
        for t, length in filter_soft_clipped(read.cigartuples):
            if t == MATCH:
                for i in xrange(length):
                    if reference_sequence[reference_index + i] != read_sequence[read_index + i] and reference_sequence[
                       reference_index + i] != "N" and read_sequence[read_index + i] != "N":
                        mismatches[reference_sequence[reference_index + i]][read_sequence[read_index + i]] += 1
                reference_index += length
                read_index += length
            if t == INSERTION:
                for i in xrange(length):
                    mismatches[None][read_sequence[read_index + i + 1]] += 1
                read_index += length
            if t == DELETION:
                for i in xrange(length):
                    mismatches[reference_sequence[reference_index + i + 1]][None] += 1
                reference_index += length

    print("*\t{0}\t{1}\t{2}\t{3}".
          format(mismatches["A"]["C"], mismatches["A"]["T"], mismatches["A"]["G"], mismatches["A"][None]))
    print("{0}\t*\t{1}\t{2}\t{3}".
          format(mismatches["C"]["A"], mismatches["C"]["T"], mismatches["C"]["G"], mismatches["C"][None]))
    print("{0}\t{1}\t*\t{2}\t{3}".
          format(mismatches["G"]["A"], mismatches["G"]["C"], mismatches["G"]["T"], mismatches["G"][None]))
    print("{0}\t{1}\t{2}\t*\t{3}".
          format(mismatches["T"]["A"], mismatches["T"]["C"], mismatches["T"]["G"], mismatches["T"][None]))
    print("{0}\t{1}\t{2}\t{3}\t*".
          format(mismatches[None]["A"], mismatches[None]["C"], mismatches[None]["T"], mismatches[None]["G"]))


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
