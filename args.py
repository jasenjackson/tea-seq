import argparse
import os
import sys

from io_utils import list_dir_complete, read_params_file, if_not_dir_make, read_feature_map
from io_utils import retrieve_paired_end_reads
from feature_search import kmers_k


def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', help='Directory containing paired end reads')
    parser.add_argument('-flash', default='flash',
                        help='Path to flash program')
    parser.add_argument(
        '-w', default=8, help='Set word size for kmer searches')
    parser.add_argument(
        '-r', help='Path to results directory if does not exist will be made')
    parser.add_argument('-n', help='Run name')
    parser.add_argument(
        '-f', help='Path to (currently) python file or dir containing list of FEATURES')
    parser.add_argument('--filter', default=False,
                       help='If the directory given in -d has multible paired \
                       end reads but you wish to use only one set pass the \
                       indentifying prefix of that set here')
    parser.add_argument('-m', help='If a directory is given to -f need a csv file \
                        that maps the indentifier prefix to a specific feature \
                        (params) file.')
    args = parser.parse_args()

    # check results dir exists and if not make

    return validate_and_process_args(args)


def validate_and_process_args(args):
    args = check_features(args)
    validate_results_paths(args)
    args = process_features(args)
    

    return args

def check_features(args):
    if os.path.isdir(args.f):
        args.m = read_feature_map(args.m)
    return args

def apply_logic(args):
    '''
    Controls which args can be given with which and makes sure user gives
    enough info to run the program.
    
    EH
    '''
    pass
    


def process_features(args):
    '''
    Reads in features from the given csv file.
    processes kmers using kmers_k method.
    '''
    
    args.f = read_params_file(args.f)
    for seq in args.f:
        # add kmers and k to sequence
        sequence = seq[1]
        kmers = kmers_k(sequence, args.w)
        seq.append(kmers)

    return args


def validate_results_paths(args):
    '''
    Makes sure that the results dir is ready to go. If it does not exist
    creates it along the run dir that is a subdir of the results dir.
    Also gets all fastq files in the data directory given by -d.
    '''
    if not os.path.exists(args.r):
        os.makedirs(args.r)
    run_dir = os.path.join(args.r, args.n)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    if args.d:
        args.d = retrieve_paired_end_reads(args.d)
            
            
