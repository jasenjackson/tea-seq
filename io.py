#!/usr/bin/python
'''
Title: file input/output functions (io.py)
Author: Jasen M. Jackson, Loyola '19
This script contains definitions that handle file i/o
'''
import os

ALLOWED_FASTA = {'.fa', '.fasta', '.fastq'}


def type_fastq(path):
    if path[-6:].lower() == '.fastq':
        return True
    else:
        return False


def fasta_check(path):
    ext = os.path.basename(path).split('.')[-1]
    # get just basename of file then split at . to create two strings
    # then just take the string after the . in ex test.fasta to
    # return 'fasta'
    if ext in ALLOWED_FASTA:
        return True
    else:
        return False
