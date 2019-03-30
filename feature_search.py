#!/usr/bin/python
'''
Title: Main pipeline wrapper for TEA-seq (main.py)
Author: Jasen M. Jackson, Loyola '19
Date: 2/20/19-
This script contains functions for sequence feature search algorithms used by TEA seq.
'''

#Return a list of Kmers that cover "threshold" % of the primerss
def kmers_threshold(primer, threshold):

	k = int(round(threshold*len(primer))) #k_min = 0.9*primer.length
	kmers = []
	for i in range(len(primer)-k+1):
		kmer = primer[i:i+k]
		kmers.append(kmer)
	return(kmers, k)

#Returns a list of Kmers of length k
def kmers_k(primer, k):
	k = int(k)
	kmers = []
	for i in range(len(primer)-k+1):
		kmer = primer[i:i+k]
		kmers.append(kmer)
	return(kmers)

#Returns the hamming distance between two strings
def hamming(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
