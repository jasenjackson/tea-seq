#!/usr/bin/python
'''
Title: Library filter functions (library_filter.py)
Author: Jasen M. Jackson, Loyola '19
This script contains functions necessary for filtering TRAPseq (teaseq???) libraries.
'''
import glob, os
import os.path
import sys
from feature_search import *
from collections import defaultdict

def merge_reads(r1,  r2, library_name, run_name, bin, data):

	#If contig has already been merged, skip this step.
	if os.path.exists('results/'+run_name+"/"+library_name+'/out.extendedFrags.fastq'):
		print("\t"+library_name+' has already been merged. Skipping to next step!')

	#otherwise merge forward and reverse read files for the library
	else:
		#run FLASH from the bin directory with 250bp maximum overlap AND interleaved output
		output_directory = "results/"+run_name+"/"+library_name
		print("\tMerging forward and reverse reads for " + library_name+"...")
		os.system(bin + "/./flash " +data+"/"+r1+ " " +data+"/"+r2
		+ " -d " + output_directory + " -M 250 --interleaved-output > results/"+
		run_name+"/"+library_name+"/flash.log")
		print("\t\tPair-end alignment console output directed to results/"+run_name+"/"+library_name+"/flash.log")

def collate(library_name, run_name):

	#If combined file has already been made, skip this step
	if os.path.exists('results/'+run_name+'/'+library_name+'/'+library_name+".combined.fastq"):
		print("\t"+library_name+' already has a combined fastq file. Skipping to next step!')

	else:
		print("\tCollating assembled and unassembled reads into combined file..")
		#create file that combines assembled and unassembled reads
		combined_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".combined.fastq"
		combined_file = open(combined_file_path, 'a')

		#dump all of the contents of the assembled file into the combined file
		assembled_file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
		assembled_file = open(assembled_file_path, 'r')
		for line in assembled_file.readlines():
			combined_file.write(line)
		assembled_file.close()
		print("\t\t"+assembled_file_path+" was succesfully added to "+combined_file_path)

		#dump all of the contents of unassembled into combined files
		unassembled_file_path = "results/"+run_name+'/'+library_name+"/out.notCombined.fastq"
		unassembled_file = open(unassembled_file_path, 'r')
		for line in unassembled_file.readlines():
			combined_file.write(line)
			unassembled_file.close()
		print("\t\t"+unassembled_file_path+" was succesfully added to "+combined_file_path)

def feature_trim(features, library_name, run_name): #feature count, filture, trim

	#open library FASTQ file & create trimmed file
	file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
	library_file = open(file_path, 'r')
	print("\tTrimming "+library_name+" with provided feature sequences...")
	trimmed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.fastq"
	if os.path.exists(trimmed_file_path): return 0
	else:
		trimmed_file = open(trimmed_file_path, 'a')

		#iterate through every line in the library sequence file and compare to every filter sequence
		has_adapter, has_element, has_killSequence, is_trimmedLine = False, False, False, False
		trimmed_line = ""
		count, adapter_pos, adapter_dist, element_pos, element_dist = 0, -1, -1, -1, -1
		while True: #for each line in the file
			header = library_file.readline()
			sequenceLine = library_file.readline()
			plusLine = library_file.readline()
			qualityScore = library_file.readline()
			if not qualityScore: break
			has_adapter, has_element, has_killSequence = False, False, False

			#interrogate sequenceLine for important sequences
			for feat in features:
				feature, threshold, type = feat[1], feat[2], feat[3]
				if (type=="remove") and (feature in sequenceLine):
					has_killSequence = True
					break
				if (type=="adapter"): # what ab RC adapter? Should be reverse orientation..
					adapter_pos, adapter_dist = slidingSearch(sequenceLine, feature, 6, 2)
					adapter_len = len(feature)
					if (adapter_pos != -1): has_adapter = True
				if (type=="element"):
					element_pos, element_dist = slidingSearch(sequenceLine, feature, 3, 1)
					if (element_pos != -1):
						has_element = True

			# trim and store eligible sequences
			if ((has_killSequence==False) and (has_adapter==True) and (has_element==True)):
				count += 1
				#print(count)
				adapter_end = adapter_pos+adapter_len
				trimmed_line = sequenceLine[adapter_end:element_pos] + '\n'
				#qualityScore = qualityScore[adapter_end:element_pos] + '\n'
				fasta_header = ">"+header
				new_entry = fasta_header+trimmed_line
				trimmed_file.write(new_entry)

		#print console output for trimmed reads process
		features_used = ""
		for feat in features:
			type = feat[3]
			if (type=="remove" or type == "adapter" or type == "element"):
				features_used = features_used + feat[0] + ", "
		features_used = features_used[:-2]
		print("\t\t"+file_path+" was succesfully trimmed and filtered using: "+ features_used)
		print("\t\t"+ str(count) +" trimmed reads added to " +trimmed_file_path)

def remove_duplicates(library_name, run_name):
	## requires
	## open trimmed file & remove duplicates
	trimmed_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".trimmed.fastq"
	trimmed_file = open(trimmed_file_path, 'r')
	print("\tRemoving dups from "+library_name+".trimmed.fastq...")
	duplicates_removed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.duplicates_removed.fastq"
	duplicates_removed = open(duplicates_removed_file_path, 'a')

	## stores non-redundant integrations in nested dictionary
	nested_dict = lambda: collections.defaultdict(nested_dict)
	redundancy_map = nested_dict()

	while True:
		header = trimmed_file.readline()
		seq = trimmed_file.readline()
		if not seq: break
		#add_to_map(seq, redundancy_map)
