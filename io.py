#!/usr/bin/python
'''
Title: file input/output functions (io.py)
Author: Jasen M. Jackson, Loyola '19
This script contains definitions that handle file i/o
'''
import os.path

def type_fa(path):
    ## boolean: path suffix is '.fasta'
    if path[-3:].lower() == ".fa":
        return True
    else: return False

def type_fasta(path):
    if path[-6:].lower() == ".fasta":
        return True
    else: return False

def type_fastq(path):
    if path[-6:].lower() == '.fastq':
        return True
    else: return False
def fasta_check(path):
    ## check file for .fasta/.fa ending, get name of genome
    fasta_check = type_fasta(path)
    fa_check = type_fa(path)
    if (fasta_check or fa_check):
        return True
    else: return False
