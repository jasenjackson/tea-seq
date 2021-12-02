# holds functions for blasting results against reference genome
import math
import subprocess
import os
import csv

from Bio import SeqIO

# could be a good idea to have a class for this but what about having
# a class that holds the genes for earier lookup?
# for right now just have gene search return two lists, first being the
# row of the blast search and the second being the row of the gene the blast
# hit is adjacent to

# TODO Working everything into a class, possibly
# could mean want to have a way to hold genes in memory in a class
# would need a way to read them in in GeneJunction

# RESULTS IN BLAST FILE CONTAINED EMPTY ENTRIES HEADER NO SEQ

BDB = '/home/ethan/Documents/Glycine_Max_RefSeq/genome_assemblies/GM_BDB/GM_2.0_BDB'
GENES = '/home/ethan/Documents/GM_genes/gene_result.txt'
OUTPUT = '/home/ethan/Documents/Illumina_files/blast_test'
FASTA = '/home/ethan/Documents/Illumina_files/results/test_2/HL-blast sequence /test_2_trimmed.duplicates_removed.fasta'
acc = '/home/ethan/Documents/Illumina_files/acc_2_chr.txt'
BLAST_results = '/home/ethan/Documents/Illumina_files/blast_test/BLAST_results'


def run_blast(output_dir, filtered_fastq, BDB, dir_name='BLAST_results'):
    '''
    Given a fasta formated file of flanking sequences identified by methods
    in the library builder file BLASTs these sequences against a local database
    and returns the path to the results. Results are output format 10 or tab
    delim.


    Gene download from NCBI header
    tax_id	Org_name	GeneID	CurrentID	Status	Symbol
    Aliases	description	other_designations	map_location
    chromosome	genomic_nucleotide_accession.version
    start_position_on_the_genomic_accession
    end_position_on_the_genomic_accession
    orientation	exon_count
    OMIM

    '''
    # need to figure out what format the query file will be coming in
    results_path = os.path.join(output_dir, dir_name)
    cmd = ['blastn', '-db', BDB, '-query',
           filtered_fastq, '-outfmt', '10', '-out', results_path]
    print(' '.join(cmd))
    subprocess.call(cmd)

    return results_path


def gene_search(genes_file, blast_results, chr_2_acc_path, gene_chr_index=10, gene_bp_index=12,
                alignment_chr_index=2, alignment_bp_index=8, max_dist=10000):
    '''
    Given a list of genes and the results from the blast search run
    in the run_blast function compares significant blast hits to the locations
    of the genes included in the genes file. The genes file should be a tab
    deliminated text file with NCBI standard formating (What you get by hitting
    the download button). 

    Returns COMING SOON
    '''
    sorted_alignments = alignment_sort(blast_results, chr_2_acc_path)
    acc_dict = acc_to_chr_dict(chr_2_acc_path)
    adjacent_pairs = []  # holds genes and seqs that were closer than threshold

    with open(genes_file) as genes:
        reader = csv.reader(genes, delimiter='\t')
        next(reader)  # skip header
        for row in reader:
            try:
                gene_chr, gene_bp = int(
                    row[gene_chr_index]), int(row[gene_bp_index])
            except ValueError:
                continue  # some genes on scaffolds do not have chr number
            search_result = nearest_binary_search(gene_chr, gene_bp,
                                                  sorted_alignments)
            if search_result:
                try:
                    aligned_seq = sorted_alignments[gene_chr-1][search_result]
                    adjacent_pairs.append((row, aligned_seq))
                except IndexError:
                    continue
                # exception for some weird cases probably write to somewhere

    return adjacent_pairs

    # do some output formating
    # else no hit continue on to the next gene

    # sort genes by position and sort hits by position then compare them
    # would allow for binary type searches within each chrosome
    # LOVE THE LOG BABY!!!!!!


def nearest_binary_search(query_chr, query_bp, sorted_genome_list, bp_index=8,
                          max_dist=10000):
    '''
    Takes in two integers, the chromsome number of the query gene and its
    base pair position. The third parameter is the sorted list of allignments
    coming from genome_sort function which ultimately gets its input material
    from BLAST allignment results. This function uses a binary search approach
    to find if there are any allignments within 10kbp (max_dis default) of the
    given query.
    '''
    bp_index = int(bp_index)
    search_space = sorted_genome_list[query_chr-1]  # chr given in base 1
    low, high, midpoint = 0, len(search_space)-1, None

    while low <= high:
            # print(low, high)
        midpoint = math.floor((low + high) / 2)
        # print(midpoint)
        # use floor to round midpoints of ls with even num of values
        midpoint_value = int(search_space[midpoint][bp_index])
        if query_bp < midpoint_value:
            high = midpoint - 1
        elif query_bp > midpoint_value:
            low = midpoint+1
        else:
            break

    closest_value = int(search_space[midpoint][bp_index])
    if abs(query_bp - closest_value) <= max_dist:
        return midpoint
    else:
        return False

# Need to convert accessions to chromosome numbers
# main problem seems to be with scaffolds where chromsome number
# is not a thing test on HL12 because it was having a lot of issues

def alignment_sort(file_to_be_sorted, chr_2_acc_path, acc_index=1, bp_index=8,
                   delim=',', header=False):
    '''
    Given a delimited file (default = \t) and the index in a row (base 0) of
    where the chromosome number and bp position within that chromosome can
    be found returns a list of lists where the index of the outer list + 1
    corresponds to the chromosome number and the inner list contains the rows
    of the original file sorted by base pair position.

    As an example writing the_list[0][-1] would return a list corresponding to
    the last element on chromosome 1. This allows for locating the nearby genes
    faster because if we know the chromsome a gene is on we can use a binary
    search approach.
    '''
    acc_dict = acc_to_chr_dict(chr_2_acc_path)

    sorted_dict = {}  # intial data structure is dictionary because we might
    # not know how many chromosomes there are. After everything is added we can
    # convert to a list for better storage without running the risk of index
    # errors

    with open(file_to_be_sorted) as ftbs:
        reader = csv.reader(ftbs, delimiter=delim)
        if header:
            next(reader)  # skip the header if present
        for row in reader:
            try:
                chr = acc_dict[row[acc_index]]
            except KeyError:
                continue
            if chr in sorted_dict:
                sorted_dict[chr].append(row)
            else:
                sorted_dict[chr] = [row]
        # all elements added under a key that == their chromsome number
        # now need to sort the elements in each chromosome

    sorted_list = [[]] * len(sorted_dict)
    for chr in sorted_dict:
        try:
            sorted_list[chr-1] = sorted(
                sorted_dict[chr], key=lambda x: int(x[bp_index]))
        except IndexError:
            print('EXCEPTION')
            for chromo, values in sorted_dict.items():
                print(print(chromo), len(values))
                print('=======================')
                print('BP Index: {} Chromosome: {}'.format(bp_index, chr))
            continue

    return sorted_list


def acc_to_chr_dict(acc_to_chr_path, delim='\t', header=True):
    with open(acc_to_chr_path) as acc:
        reader = csv.reader(acc, delimiter=delim)
        if header:
            next(reader)
        return {acc: int(chromo) for chromo, acc in reader}

#s = alignment_sort(BLAST_results, acc)
