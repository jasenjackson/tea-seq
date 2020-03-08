#!/usr/bin/python
'''
Title: Deprecated definitions (library_filter.py)
Author: Jasen M. Jackson, Loyola '19
This script contains deprecated functions from TRAPseq (teaseq???) workflow.
'''
def feature_count(features, library_name, run_name):

	#open library FASTQ file
	file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
	library_file = open(file_path, 'r')
	print("\tSearching "+library_name+" for provided filter sequences...")

	#create feature_count file
	feature_count_file_path = "results/"+run_name+"/"+library_name+"/feature_count.csv"
	feature_count_file = open(feature_count_file_path, 'a')
	feature_table_entry=""

	#iterate through each filter sequence and add all matching reads to new file
	for seq in features:

		#get name, sequence, match threshold, kmers and appropriate k parameter from each feature
		feature_name, sequence, match_threshold = seq[0], seq[1], seq[2]
		kmers, k = kmers_threshold(sequence, match_threshold)

		#iterate through every line in the library sequence file
		i = 0 #line number
		feature_count = 0
		for line in library_file:
			i = i+1
			if ((i+2)%4 == 0): #i.e: for every "Sequence" line
				for kmer in kmers: #search all kmers
					if kmer in line: #if kmer is in the read
						feature_count = feature_count + 1
						break
		library_file.seek(0) #return pointer to the top of the library file, count next feature

		#compute percentage of reads that matched
		total_reads = (i/4)
		if total_reads > 0: #avoid dividing by zero
			feature_percentage = '{:0.3f}'.format((float(feature_count)/float(total_reads)))
		else:
			feature_percentage = "0%"

		#Primer: X reads (0%)
		print("\t\t"+feature_name+": "+str(feature_count)+" reads ("+str(feature_percentage)+") k="+str(k))

		#add to row to feature table csv
		feature_table_entry = feature_table_entry + str(feature_count) + "," + str(feature_percentage) + ","

	#add entire feature row for library
	print("\t\tFeature table entry: " + feature_table_entry)
	feature_count_file.write(feature_table_entry+'\n')

def remove_duplicates2(library_name, run_name):
		#open trimmed file & remove duplicates
		trimmed_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".trimmed.fastq"
		trimmed_file = open(trimmed_file_path, 'r')
		print("\tRemoving dups from "+library_name+".trimmed.fastq...")
		duplicates_removed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.duplicates_removed.fastq"
		duplicates_removed = open(duplicates_removed_file_path, 'a')

		unique_set = [] #collect unique seqs
		duplicate_count = 0

		while True:
			header = trimmed_file.readline()
			sequenceLine = trimmed_file.readline()
			if not sequenceLine: break

			# compare last 20bp of sequences to "seen" set
			sequence_end = sequenceLine[-20:] # make variable

			if sequence_end not in unique_set: # TODO: compare hamming distances
				duplicates_removed.write(header+sequenceLine)
				unique_set.append(sequence_end)
			else:
				duplicate_count += 1
		print(unique_set[1:5])
		print("\t\tdups: " + str(duplicate_count))
		print("\t\tunique: " + str(len(unique_set)))
