__author__ = 'nikita_kartashov'

from sys import argv
from os import path, remove, rename, makedirs
from Bio import SeqIO
from Bio.SeqUtils import GC

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def relative_path(relative_path_part):
    script_dir = path.dirname(__file__)
    return path.join(script_dir, relative_path_part)

DEFAULT_INPUT = list(map(relative_path, ['../../resource/hw1/test.fastq']))

RESULT_PATH = relative_path('../../results/hw1')


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


def get_record_quality(record):
    return record.letter_annotations["phred_quality"]


def get_gc_content(records):
    # Precision is 0.1
    BINS = 100 * 10
    gc_contents = [GC(r.seq) for r in records if not is_read_low_quality(get_record_quality(r))]
    plt.hist(gc_contents, BINS)
    plt.ylabel('Number of reads')
    plt.xlabel('%GC')
    plt.savefig(result_file('gc_contents.png'))
    plt.clf()


def combine_two_scores(old, new):
    if len(old) < len(new):
        old, new = new, old
    for i, e in enumerate(new):
        old[i] = (old[i] + e) * 1.0 / 2
    return old


def score_to_probabilites(scores):
    return list(10 ** ((-score) * 1.0 / 10 + 2) for score in scores)


def get_quality_distribution(records):
    qualities = map(get_record_quality, records)
    mean_scores = reduce(combine_two_scores, qualities, [])
    probabilities = score_to_probabilites(mean_scores)
    plt.plot(range(0, len(probabilities)), probabilities)
    plt.ylabel('Error probability %')
    plt.xlabel('Base')
    plt.savefig(result_file('nucleotide_error_probability.png'))
    plt.clf()

if __name__ == '__main__':
    make_result_dir()

    argv = argv[1:]
    if len(argv) == 0:
        argv.extend(DEFAULT_INPUT)

    for filename in argv:
        # preprocess_file(filename)
        with open(filename) as input_file:
            records = get_records(input_file)
            get_gc_content(records)
            records = get_records(input_file)
            get_quality_distribution(records)







