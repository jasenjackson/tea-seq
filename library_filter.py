#!/usr/bin/python
'''
Title: Library filter functions (library_filter.py)
Author: Jasen M. Jackson, Loyola '19
This script contains functions necessary for filtering TRAPseq (teaseq???) libraries.
'''
import glob
import os
import os.path
import sys
from io import *
from redundancy_map import *
from feature_search import *
from collections import defaultdict

from io import if_not_dir_make


def feature_trim(features, library_name, run_name, end_size):  # feature count, filture, trim

    # open library FASTQ file & create trimmed file
    file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"  # extended file
    library_file = open(file_path, 'r')
    print("\tTrimming "+library_name+" with provided feature sequences...")
    trimmed_file_path = "results/"+run_name+"/" + \
        library_name+"/"+library_name+".trimmed.fasta"
    if os.path.exists(trimmed_file_path):
        return 0
    else:
        trimmed_file = open(trimmed_file_path, 'a')

        # iterate through every line in the library sequence file and compare to every filter sequence
        has_adapter, has_element, has_killSequence, is_trimmedLine = False, False, False, False
        trimmed_line = ""
        count, adapter_pos, adapter_dist, element_pos, element_dist = 0, -1, -1, -1, -1
        while True:  # for each line in the file
            header = library_file.readline()
            sequenceLine = library_file.readline()
            plusLine = library_file.readline()
            qualityScore = library_file.readline()
            if not qualityScore:
                break
            has_adapter, has_element, has_killSequence = False, False, False

            # interrogate sequenceLine for important sequences
            for feat in features:
                feature, threshold, type = feat[1], feat[2], feat[3]
                if (type == "remove") and (feature in sequenceLine):
                    has_killSequence = True
                    break
                # what ab RC adapter? Should be reverse orientation..
                if (type == "adapter"):
                    adapter_pos, adapter_dist = slidingSearch(
                        sequenceLine, feature, 6, 2)
                    adapter_len = len(feature)
                    if (adapter_pos != -1):
                        has_adapter = True
                if (type == "element"):
                    element_pos, element_dist = slidingSearch(
                        sequenceLine, feature, 3, 1)
                    if (element_pos != -1):
                        has_element = True

            # trim and store eligible sequences
            if ((has_killSequence == False) and (has_adapter == True) and (has_element == True)):
                # print(count)
                adapter_end = adapter_pos+adapter_len
                trimmed_line = sequenceLine[adapter_end:element_pos] + '\n'
                if len(trimmed_line) >= end_size:
                    count += 1
                    fasta_header = ">"+header
                    new_entry = fasta_header+trimmed_line
                    trimmed_file.write(new_entry)

        # print console output for trimmed reads process
        features_used = ""
        for feat in features:
            type = feat[3]
            if (type == "remove" or type == "adapter" or type == "element"):
                features_used = features_used + feat[0] + ", "
        features_used = features_used[:-2]
        print("\t\t"+file_path +
              " was succesfully trimmed and filtered using: " + features_used)
        print("\t\t" + str(count) + " trimmed reads added to " + trimmed_file_path)


def remove_duplicates(library_name, run_name, end_size):

    # create redundancy map
    trimmed_file_path = "results/"+run_name+'/' + \
        library_name+"/"+library_name+".trimmed.fasta"
    rm = redundancy_map(trimmed_file_path, end_size)
    print("\t\tUnique reads: "+str(rm.unique_count))
    print("\t\tDuplicate reads: "+str(rm.duplicate_count))
    print("\t\tTotal reads: " + str(rm.total_count))
    rm.print_map_head(10)
    '''## input/output function
	## open trimmed file & remove duplicates
	trimmed_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".trimmed.fastq"

	#duplicates_removed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.duplicates_removed.fastq"
	#duplicates_removed = open(duplicates_removed_file_path, 'a')
	print("\tRemoving dups from "+library_name+".trimmed.fastq...")

	## stores non-redundant integrations in nested dictionary
	nested_dict = lambda: defaultdict(nested_dict)
	redundancy_map = nested_dict()

	while True:
		header = trimmed_file.readline()
		seq = trimmed_file.readline()
		if not seq: break
		#add_to_map(seq, redundancy_map)
	'''
