#!/usr/bin/python
'''
Title: Main pipeline wrapper for TEA-seq (main.py)
Author: Jasen M. Jackson, Loyola '19
Date: 2/20/19-
This script contains the main workflow for the TEA-seq pipeline.
'''
from args import set_args
from library_builder import make_libraries

def main():
    argue = set_args()
    make_libraries(argue.r, argue.n, argue.d, argue.f, argue.flash)
    

if __name__ == "__main__":
    main()

