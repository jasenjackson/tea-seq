import argparse
import os

from io import list_dir_complete, read_params_file, if_not_dir_make
from feature_search import kmers_k


def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', help='Directory containing paired end reads')
    parser.add_argument('-t', default='./bin',
                        help='Path to temp dir for short term file storage')
    parser.add_argument(
        '-w', default=8, help='Set word size for kmer searches')
    parser.add_argument(
        '-r', help='Path to results directory if does not exist will be made')
    parser.add_argument('-n', help='Run name')
    parser.add_argument(
        '-f', help='Path to (currently) python file containing list of FEATURES')
    parser.add_argument('--filter', default=False,
                       help='If the directory given in -d has multible paired \
                       end reads but you wish to use only one set pass the \
                       indentifying prefix of that set here')
    args = parser.parse_args()

    # check results dir exists and if not make

    return validate_and_process_args(args)


def validate_and_process_args(args):
    validate_paths(args)
    args = process_features(args)

    return args


def process_features(args):
    args.f = read_params_file(args.f)
    for seq in args.f:
        # add kmers and k to sequence
        sequence = seq[1]
        kmers = kmers_k(sequence, args.w)
        seq.append(kmers)

    return args


def validate_results_paths(args):
    if not os.path.exists(args.r):
        os.makedirs(args.r)
    run_dir = os.path.join(args.r, args.n)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    if args.d:
        args.d = list_dir_complete(args.d)


def make_libraries(args):
    run_path = os.path.join(args.r, args.n)
    for identifier, fastq_files in args.d.items:
        library_path, made_dir = if_not_dir_make(run_path, identifier)
