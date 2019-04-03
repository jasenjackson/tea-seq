#!/usr/bin/python
'''
Title: Main pipeline wrapper for TEA-seq (feature_search.py)
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

def kmer_search_depr(query, kmers, subject, threshold, direction):
	min_length = threshold*len(query)
	pos, offset, dist = 0, 0, 0
	match, subQuery = "", ""

	#search for location beginning w/ first kmer in kmers list
	for kmer in kmers:
		pos = subject.find(kmer)
		if pos != -1 and direction=='f':
			subQuery = query[offset:]
			match = subject[pos:pos+len(subQuery)]
			dist = hamming(match, subQuery)
			if (dist <= MAXDIST and len(match) >= min_length):
				return pos, pos+len(subQuery)
		if pos != -1 and direction=='r':
			kmers.reverse() # assumes left-to-right kmer orientation
			subQuery = query[:len(query)-offset]
			match_start = pos+len(kmer)-len(subQuery)
			match_end = pos+len(kmer)
			match = subject[match_start:match_end]
			dist = hamming(match, subQuery)
			if (dist <= MAXDIST and len(match) >= min_length):
				return match_start, match_end
		offset += 1
	# if no qualifying match is found
	return -1, -1


#Returns starting position of best subject-query alignment.
def kmer_search(subject, query, maxdist, breakdist):
    k = len(query) # kmers should be query length
    curr_dist, match_start, lowest_obs_dist = 0, 0, k
    for i in range(0, len(subject)-k):
        curr_dist = hamming(subject[i:i+k], query)
        if curr_dist < lowest_obs_dist: # update best match
            lowest_obs_dist = curr_dist
            match_start = i
            if curr_dist <= breakdist: # Heuristic: kill for exact match
                break
    if lowest_obs_dist <= maxdist: return(match_start, lowest_obs_dist)
    else: return(k, k)
