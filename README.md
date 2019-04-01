# Targeted Sequencing of Retrotransposon Integrations in Soybean populations  [work in progress]
This document contains the documentation for the analysis of transposable-element anchored PCR sequencing (tea-seq) data obtained by the method described in the following document: https://bit.ly/2WzWjmg. 
<center><a href = "https://github.com/jasenjackson/tea-seq/"><img src = "https://github.com/jasenjackson/tea-seq/blob/master/Screen%20Shot%202019-03-31%20at%204.18.19%20PM.png?raw=true"/></a></center>

Repetitive elements predominate the landscape of virtually all plant & animal genomes. Retrotransposons are especially ubiquitous because of their ability to rapidly amplify their copy-number in eukaryotic genomes via reverse-transcription of their RNA intermediate. Because of their relative abundance and ability to regulate the expression of nearby genes, it has been proposed that some of these elements have been co-opted into mutualistic relationships by their hosts, whereby, under immense selective pressure, transcriptionally-silenced retroelements are unleashed into coding regions, resulting in a rapid diversification of the host’s genome. Genome-wide analysis of retrotransposon insertional polymorphisms between wild and cultivated soybean may reveal novel retroelement insertions in the regulatory regions of domestication-related traits. As such, we have developed a transposon-anchored PCR protocol which amplifies genomic regions flanking target retrotransposon families and subjects their amplicon libraries to Illumina sequencing. 

## Overview
Our computational strategy depends on cleaning and annotating paire

0. modify rm_feats to return list of features [done]
1. Convert rm into class [done]
2. Count duplicates & uniques, store as rm attribute [done]
3. Use dups on trimmed_file [done]
4. Align trimmed_dups_removed file to reference without LTR added [in progress]
5. Fix quality/adapter trimming [in progress]
    a. quality TRIM paired-end reads (Trimmomatic), assess quality graphs before & after (FASTQC)
    b. merge reads before and after quality trimming, assess quality (FASTQC)
    c. trim adapter from 5'/3' end (Trimmomatic)
6. Clean up genome-junction regions w/ redundancy map [in progress]
    d. locate LTR edge using KMER search. Create duplicate file with LTRs removed.
    e. create redundancy map from trimmed file (extend to paired-end read), remove polypurine tract sequences
6. Align *with* LTR added, allow for no more than <3> multi-hits. Remove solid hits.
7. Align *w/o* LTR, allow for no more than <10> hits. Remove repeats. Characterize
8. Characterize locations by reading BAM file. (look in GENEIOUS, use pysam to fetch regions)
9. Compare to known locations.
10. Repeat for all chromosomes & high quality 

## Cleaning up the FASTQ libraries
The first step in the analysis is to process the raw paired-end read FASTQ files for  for each amplicon library. This can be accomplished using the 'library_filter.py' script. This is still a work in progress, but it has multiple steps:
1. Merge the paired-end reads using FLASH.
2. Collect merged and un-merged reads into a single file (incomplete)
3. Trim sequences between the adapter and the end of the LTR.
4. Filter out duplicates and sequences containing the poly-purine tract

The 'library_filter' script takes one argument: run_name. Each time the program is executed it creates a different folder in the results file, using the run_name variable.
        
        python library_filter.py "run name"

## Identify flanking sequences in other repetitive elements
Next, we blast each amplicon library to the soybean transposable element database (https://www.soybase.org/soytedb/) and generate a list of non-redundant, non-nested flanking sequences. To do this, you can download the reference genome as a FASTA file from the SoyTEDB website and run the following command to convert the raw fasta file into a BLAST database. 

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

![non-nested gMR30](https://github.com/jasenjackson/tea-seq/blob/master/non-nested-GMR30-sites-HL2.png?raw=true)
      
## Generate a list of "known" locations for your transposon family of interest
It would be helpful to know if our amplicon libraries contain all of the previously annotated family members. It is also very likely that the current reference genome does not have every family member, and we would like to know if our methods have uncovered new members/integration sites. We have a database containing the names and positions of each member of each retrotransposon family (with sequences and descriptions), but the positions are for the old soybean reference assembly. We need to translate these positions to the new assembly! 




