#!/usr/bin/python
'''
Title: Main pipeline wrapper for TEA-seq (main.py)
Author: Jasen M. Jackson, Loyola '19
Date: 2/20/19-
This script contains the main workflow for the TEA-seq pipeline.
'''
from args import set_args
from library_builder import make_libraries

from blast import run_blast, gene_search
from io_utils import write_adjacent_pairs

def main():
    argue = set_args()
    filtered_fastas = make_libraries(
        argue.r, argue.n, argue.d, argue.f, argue.flash)
    # now need to do the blasting for each plant's filtered fasta
    
    for library_path, filtered_fasta in filtered_fastas:
        print(library_path, filtered_fasta)
        alignment_path = run_blast(library_path, filtered_fasta, argue.b)
        adjacent_pairs = gene_search(argue.g, alignment_path, argue.a)
        write_adjacent_pairs(adjacent_pairs, library_path)
    


if __name__ == "__main__":
    main()
