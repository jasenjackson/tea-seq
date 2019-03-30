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

def create_paths(run_name):

	# create results directory if it does not exist.
	if not os.path.exists("results"):
		os.makedirs("results")

	# create results/run_name/ directory, if it does not exist
	if not os.path.exists('results/'+run_name):
		os.makedirs('results/'+run_name)

def merge_reads(r1,  r2, library_name, run_name):

	#If contig has already been merged, skip this step.
	if os.path.exists('results/'+run_name+"/"+library_name+'/out.extendedFrags.fastq'):
		print(lt+"\t"+library_name+' has already been merged. Skipping to next step!')

	#otherwise merge forward and reverse read files for the library
	else:
		#run FLASH from the bin directory with 250bp maximum overlap AND interleaved output
		output_directory = "results/"+run_name+"/"+library_name
		print(lt+"\tMerging forward and reverse reads for " + library_name+"...")
		os.system(BIN_DIR + "/./flash " +DATA_DIR+"/"+r1+ " " +DATA_DIR+"/"+r2
		+ " -d " + output_directory + " -M 250 --interleaved-output > results/"+
		run_name+"/"+library_name+"/flash.log")
		print(lt+"\t\tPair-end alignment console output directed to results/"+run_name+"/"+library_name+"/flash.log")

def collate(library_name, run_name):

	#If combined file has already been made, skip this step
	if os.path.exists('results/'+run_name+'/'+library_name+'/'+library_name+".combined.fastq"):
		print(lt+"\t"+library_name+' already has a combined fastq file. Skipping to next step!')

	else:
		print(lt+"\tCollating assembled and unassembled reads into combined file..")
		#create file that combines assembled and unassembled reads
		combined_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".combined.fastq"
		combined_file = open(combined_file_path, 'a')

		#dump all of the contents of the assembled file into the combined file
		assembled_file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
		assembled_file = open(assembled_file_path, 'r')
		for line in assembled_file.readlines():
			combined_file.write(line)
		assembled_file.close()
		print(lt+"\t\t"+assembled_file_path+" was succesfully added to "+combined_file_path)

		#dump all of the contents of unassembled into combined files
		unassembled_file_path = "results/"+run_name+'/'+library_name+"/out.notCombined.fastq"
		unassembled_file = open(unassembled_file_path, 'r')
		for line in unassembled_file.readlines():
			combined_file.write(line)
			unassembled_file.close()
		print(lt+"\t\t"+unassembled_file_path+" was succesfully added to "+combined_file_path)

def feature_count(features, library_name, run_name):

	#open library FASTQ file
	file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
	library_file = open(file_path, 'r')
	print(lt+"\tSearching "+library_name+" for provided filter sequences...")

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
		print(lt+"\t\t"+feature_name+": "+str(feature_count)+" reads ("+str(feature_percentage)+") k="+str(k))

		#add to row to feature table csv
		feature_table_entry = feature_table_entry + str(feature_count) + "," + str(feature_percentage) + ","

	#add entire feature row for library
	print(lt+"\t\tFeature table entry: " + feature_table_entry)
	feature_count_file.write(feature_table_entry+'\n')


def kmerSearch(query, kmers, subject, threshold, direction):
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
def slidingSearch(subject, query, maxdist, breakdist):
    k = len(query) # kmers should be query length
    curr_dist, match_start, lowest_obs_dist = 0, 0, k
    #iterate through subject with sliding window of size k=len(query)
    for i in range(0, len(subject)-k):
        curr_dist = hamming(subject[i:i+k], query) #align subject & query
        #print(subject[i:i+k]+" "+str(curr_dist))
        if curr_dist < lowest_obs_dist: # update best match
            lowest_obs_dist = curr_dist
            match_start = i
            if curr_dist <= breakdist: # Heuristic: if the match is "exceptional"
                break
    if lowest_obs_dist <= maxdist: # if the best match meets threshold requirements
        return(match_start, lowest_obs_dist)
    else:
        return(-1, -1)

def feature_trim(features, library_name, run_name): #feature count, filture, trim

	#open library FASTQ file & create trimmed file
	file_path = "results/"+run_name+'/'+library_name+"/out.extendedFrags.fastq"
	library_file = open(file_path, 'r')
	print(lt+"\tTrimming "+library_name+" with provided feature sequences...")
	trimmed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.fastq"
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
	print(lt+"\t\t"+file_path+" was succesfully trimmed and filtered using: "+ features_used)
	print(lt+"\t\t"+ str(count) +" trimmed reads added to " +trimmed_file_path)

def remove_duplicates2(library_name, run_name):
		#open trimmed file & remove duplicates
		trimmed_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".trimmed.fastq"
		trimmed_file = open(trimmed_file_path, 'r')
		print(lt+"\tRemoving dups from "+library_name+".trimmed.fastq...")
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
		print(lt+"\t\tdups: " + str(duplicate_count))
		print(lt+"\t\tunique: " + str(len(unique_set)))

def remove_duplicates(library_name, run_name):
		#open trimmed file & remove duplicates
		trimmed_file_path = "results/"+run_name+'/'+library_name+"/"+library_name+".trimmed.fastq"
		trimmed_file = open(trimmed_file_path, 'r')
		print(lt+"\tRemoving dups from "+library_name+".trimmed.fastq...")
		duplicates_removed_file_path = "results/"+run_name+"/"+library_name+"/"+library_name+".trimmed.duplicates_removed.fastq"
		duplicates_removed = open(duplicates_removed_file_path, 'a')

		unique_set = [] #collect unique seqs
		duplicate_count = 0
		redundancy_map = dict()

		while True:
			header = trimmed_file.readline()
			seq = trimmed_file.readline()
			if not seq: break

			# grab last 20bp
			new_key = seq[-20:] # make variable
			new_len = len(seq)

			# if not in redundancy map, add it
			# else, find longest
