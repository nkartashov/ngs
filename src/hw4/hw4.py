__author__ = 'nikita_kartashov'

from sys import argv
from os import makedirs, path

import pysam
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from Bio import SeqIO


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


def get_correction_rates(genome, not_corrected, corrected):
    mismatches = set()
    for read in not_corrected.fetch():
        for start, end in read.get_blocks():
            reference_sequence = genome[start: end]
            read_sequence = read.query_alignment_sequence
            for i in xrange(len(reference_sequence)):
                if read_sequence[i] != reference_sequence[i] and \
                                reference_sequence[i] != "N" and \
                                read_sequence[i] != "N":
                    mismatches.add(start + i)

    mismatches_in_corrected = set()
    for read in corrected.fetch():
        for start, end in read.get_blocks():
            reference_sequence = genome[start: end]
            read_sequence = read.query_alignment_sequence
            for i in xrange(len(reference_sequence)):
                if read_sequence[i] != reference_sequence[i] and \
                                reference_sequence[i] != "N" and \
                                read_sequence[i] != "N":
                    mismatches_in_corrected.add(start + i)

    FP = len(mismatches)
    FN = len(mismatches_in_corrected.difference(mismatches))
    TN = len(mismatches.difference(mismatches_in_corrected))
    print('FP: {0}\tFN: {1}\tTN: {2}'.format(FP, FN, TN))


if __name__ == '__main__':
    if len(argv) < 4:
        print("Reference or alignment_result was not supplied")
        exit(1)

    _, reference_path, not_corrected_path, corrected_path = argv
    reference = (list(SeqIO.parse(reference_path, 'fasta'))[0]).seq
    with pysam.AlignmentFile(not_corrected_path, "r") as not_corrected:
        with pysam.AlignmentFile(corrected_path, "r") as corrected:
            get_correction_rates(reference, not_corrected, corrected)


