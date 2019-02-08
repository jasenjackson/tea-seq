# Sequencing and Analysis of Transposable-element Insertional Polymorphisms in Soybean
This document contains the documentation for the analysis of transposable-element anchored PCR sequencing (tea-seq) data obtained by the method described in the following document:
* https://docs.google.com/document/d/1z6l_yjRtnGkOStGcOwBBOVk0gpyeZiGODGK2Vp9sF5M/edit?usp=sharing

## Cleaning up the FASTQ libraries
The first step in the analysis is to process the raw paired-end read FASTQ files for  for each amplicon library. This can be accomplished using the 'library_filter.py' script. This is still a work in progress, but it has multiple steps:
1. Merge the paired-end reads using FLASH.
2. Collect merged and un-merged reads into a single file (incomplete)
3. Trim sequences between the adapter and the end of the LTR.
4. Filter out duplicates and sequences containing the poly-purine tract

The 'library_filter' script takes one argument: run_name. Each time the program is executed it creates a different folder in the results file, using the run_name variable.
        
        python library_filter.py "run name"

## Identify flanking sequences in other repetitive elements
Next, we blast each amplicon library to the soybean transposable element database () and generate a list of non-redundant, non-nested flanking sequences. To do this, you can download the reference genome as a FASTA file from the SoyTEDB website and run the following command to convert the raw fasta file into a BLAST database. 

    cat data/soytedb.fa | \
    makeblastdb -dbtype nucl -title soytedb -out blast/soytedb

Then, to BLAST against this database, run the following command:

    cat data/soy.fa | \
    makeblastdb -dbtype nucl -title soy -out blast/soy
    ./bin/blastn \
    -query query_file.fa \
    -db blast_db_path \
    -outfmt “10 salltitles evalue” \
    -out results_path
    
The 'BLAST_parser_TE_families.py' script parses these results and stores the ID of every flanking sequence that had a significant alignment to SoyTEDB. The 'remove_sequences_from_list.py' script uses the output from 'BLAST_parser_TE_families.py' and creates an identical FASTA amplicon library with the nested flanking sequences removed. 

## Align non-nested flanking sequences to the reference genome
The genome of the Williams 82 soybean cultivar has been sequenced and uploaded to NCBI. We can align our non-nested flanking sequences directly to the reference genome. To do this, it may be easier to download the reference genome as a FASTA file, convert it into a BLAST database and use the non-nested flanking sequence FASTA file as a query. The 'BLAST_parser_ref_seq.py' script was written to parse these results. It outputs the chromosome and position of a top alignment from each non-nested flanking sequence, along with whether or not the sequence aligned within 10kbp of a gene. The 'chromosome_ideogram.R' script takes these positions and plots them, giving a visualization that looked like the following.

![alt text](https://raw.githubusercontent.com/jasenjackson/tea-seq/Screen Shot 2019-02-08 at 1.35.27 PM.png
      
## Generate a list of "known" locations for your transposon family of interest
It would be helpful to know if our amplicon libraries contain all of the previously annotated family members. It is also very likely that the current reference genome does not have every family member, and we would like to know if our methods have uncovered new members/integration sites. We have a database containing the names and positions of each member of each retrotransposon family (with sequences and descriptions), but the positions are for the old soybean reference assembly. We need to translate these positions to the new assembly! 




