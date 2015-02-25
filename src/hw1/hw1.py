__author__ = 'nikita_kartashov'

from sys import argv
from os import path, remove, rename, makedirs
from Bio import SeqIO
from Bio.SeqUtils import GC

import matplotlib.pyplot as plt


def relative_path(relative_path_part):
    script_dir = path.dirname(__file__)
    return path.join(script_dir, relative_path_part)


DEFAULT_INPUT = list(map(relative_path, ['../../resource/hw1/test3.fastq']))

RESULT_PATH = relative_path('../../resource/hw1/result')


def result_file(filename):
    return path.join(RESULT_PATH, filename)


def make_result_dir():
    if not path.exists(RESULT_PATH):
        makedirs(RESULT_PATH)


def preprocess_file(filename):
    with open(filename) as input_file:
        buffer_filename = filename + 'laofafsadf'
        with open(buffer_filename, 'w') as output_file:
            for line in input_file.readlines():
                if line.startswith('+'):
                    output_file.write('+\n')
                else:
                    output_file.write(line)
    remove(filename)
    rename(buffer_filename, filename)


def get_records(input_file):
    return SeqIO.parse(input_file, "fastq")


def is_nucleotide_low_quality(quality):
    return quality < 30


def is_read_low_quality(read_quality):
    def mapper(quality):
        return int(not is_nucleotide_low_quality(quality))

    return sum(list(map(mapper, read_quality))) * 1.0 / len(read_quality) < 0.75


def get_gc_content(records):
    # Precision is 0.1
    BINS = 100 * 10
    gc_contents = [GC(r.seq) for r in records if not is_read_low_quality(r.letter_annotations["phred_quality"])]
    plt.hist(gc_contents, BINS)
    plt.ylabel('Number of reads')
    plt.xlabel('%GC')
    # plt.show()
    plt.savefig(result_file('gc_contents.png'))


def get_quality_distribution(input):
    pass


if __name__ == '__main__':
    make_result_dir()

    argv = argv[1:]
    if len(argv) == 0:
        argv.extend(DEFAULT_INPUT)

    for filename in argv:
        preprocess_file(filename)
        with open(filename) as input_file:
            records = list(get_records(input_file))
            get_gc_content(records)







