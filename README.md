# Sequencing and Analysis of Transposable-element Insertional Polymorphisms in Soybeans
This document contains the documentation for the analysis of transposable-element anchored PCR sequencing (tea-seq) data obtained by the method described in the following document:
* https://docs.google.com/document/d/1z6l_yjRtnGkOStGcOwBBOVk0gpyeZiGODGK2Vp9sF5M/edit?usp=sharing

# Cleaning up the FASTQ libraries
The first step in the analysis is to process the raw paired-end read FASTQ files for  for each amplicon library. This can be accomplished using the 'library_filter.py' script. This is still a work in progress, but it has multiple steps:
1. Merge the paired-end reads using FLASH.
2. Collect merged and un-merged reads into a single file (incomplete)
3. Trim sequences between the adapter and the end of the LTR.
4. Filter out duplicates and sequences containing the poly-purine tract

# Identify flanking sequences in other repetitive elements
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

# Align non-nested flanking sequences to the reference genome
The genome of the Williams 82 soybean cultivar has been sequenced and uploaded to NCBI. We can align our non-nested flanking sequences directly to the reference genome. To do this, it may be easier to download the reference genome as a FASTA, convert it into a BLAST database and use the non-nested flanking sequence FASTA file as a query. 
