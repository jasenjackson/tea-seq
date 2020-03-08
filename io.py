#!/usr/bin/python
'''
Title: file input/output functions (io.py)
Author: Jasen M. Jackson, Loyola '19
This script contains definitions that handle file i/o
'''
import os
import csv
import os.path


def type_fa(path):
    # boolean: path suffix is '.fasta'
    if path[-3:].lower() == ".fa":
        return True
    else:
        return False


def type_fasta(path):
    if path[-6:].lower() == ".fasta":
        return True
    else:
        return False


def type_fastq(path):
    if path[-6:].lower() == '.fastq':
        return True
    else:
        return False


def fasta_check(path):
    # check file for .fasta/.fa ending, get name of genome
    fasta_check = type_fasta(path)
    fa_check = type_fa(path)
    if (fasta_check or fa_check):
        return True
    else:
        return False


def retrieve_paired_end_reads(read_dir):
    '''
    Returns paired end reads as a dictionary given a directory
    containing all illumina reads. Works based on the fact that
    paired reads have the same file prefix followed by _.
    Uses this prefix as the identifier and stores paired read file
    paths in a dictionary where key = indentifier and values = list
    of the two paired reads.

    EH
    '''
    paired_dict = {}
    fastq_files = [os.path.join(read_dir, fq)
                   for fq in os.listdir(read_dir) if fq[-5:] == 'fastq']
    for fq in fastq_files:
        identifier = os.path.basename(fq).split('_')[0]
        if identifier in paired_dict:
            paired_dict[identifier].append(fq)
        else:
            paired_dict[identifier] = [fq]
    return paired_dict


def if_not_dir_make(parent_dir, dir_name):
    '''
    Given a path to a parent dir and a new dir name checks if
    the complete path to the dir_name is a dir. If not makes
    that dir using os.mkdir. In either case returns the full path
    to the complete dir (parent_dir/dir_name)

    EH
    '''
    flag = False
    full_dir = os.path.join(parent_dir, dir_name)
    if not os.path.isdir(full_dir):
        os.mkdir(full_dir)
        flag = True
    return full_dir, flag

def read_params_file(path, header=True):
    '''
    Reads the params.csv file and returns it as a list of lists
    where eahc list represents one row of the params file.
    
    EH
    '''
    with open(path) as p:
        reader = csv.reader(p)
        if header:
            next(reader)
        return [row for row in reader]

def list_dir_complete(dir, filetype=None):
    if filetype:
        return [os.path.join(dir, f) for f in os.listdir(dir) if f.split('.')[-1] == filetype]
    else:
        return [os.path.join(dir, f) for f in os.listdir(dir)]
