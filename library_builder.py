import os

from io_utils import if_not_dir_make
from feature_search import kmer_search
from redundancy_map import RedundancyMap

from Bio import SeqIO
from Bio.Seq import Seq

# TODO: if each library required a different feature file will require
# a little bit of reworking. Use dictionary for lookup


def make_libraries(results_dir, run_name, fastq_dir_dict, features, flash_exe):
    run_path = os.path.join(results_dir, run_name)
    filtered_fastas = []

    for identifier, fastq_files in fastq_dir_dict.items():
        library_path, made_dir = if_not_dir_make(run_path, identifier)
        # makes new directory in run_path
        print('Merging Reads')
        combined_reads, uncombined_reads = merge_reads(
            fastq_files[0], fastq_files[1], library_path, flash_exe)

        all_reads = collate(library_path, run_name,
                            combined_reads, uncombined_reads)
        print('Trimming Reads')
        trimmed_reads = feature_trim(
            features, library_path, run_name, all_reads)
        
        filtered_fastas.append((library_path, remove_duplicates(
            library_path, run_name, trimmed_reads)))
        
    
    return filtered_fastas
        # remove duplicates returns file path to filtered fasta file
        # passes to redund map but probably will need to pull out some
        # file paths from here later


def merge_reads(r1,  r2, library_path, flash_exe):
    '''
    Given paths for pair of paired end reads fastq files (r1, r2), the path
    to a library directory and the path to the flash executable runs flash
    to merge the paired end reads into one file. Paired reads end up in
    /library_path/out.extendedFrags.fastq and uncombined reads end up in
    /library_path/out.notCombined.fastq. Returns the path to combined and
    uncombined files in a tuple in that order.

    EH
    '''
    # probably should be a way to give output paths to flash EH

    FLASH_COMBINED = os.path.join(library_path, 'out.extendedFrags.fastq')
    FLASH_UNCOMBINED = os.path.join(library_path, 'out.notCombined.fastq')

    flash_log_path = os.path.join(library_path, 'flash.log')

    if not os.path.exists(FLASH_COMBINED):
        # run flash
        cmd = '{} {} {} -d {} -M 250 --interleaved-output > {}'.format(
            flash_exe, r1, r2, library_path, flash_log_path)
        os.system(cmd)

    return FLASH_COMBINED, FLASH_UNCOMBINED


def collate(library_path, run_name, combined_flash_file, uncombined_flash_file):
    '''
    Given the library path, the name of the run and the paths to combined and
    uncombined flash result files merges the uncombined and combined results
    into a new file located at library_path/run_name_combined.fastq. Returns
    the path to the new combined file.
    '''
    combined_file_path = os.path.join(
        library_path, '{}_combined.fastq'.format(run_name))

    flash_combined_contents = open(combined_flash_file).readlines()
    flash_uncombined_contents = open(uncombined_flash_file).readlines()
    with open(combined_file_path, 'w') as cff:
        cff.writelines(flash_combined_contents)
        cff.writelines(flash_uncombined_contents)

    return combined_file_path


def feature_trim(features, library_path, run_name, extended_file, end_size=20):
    # features come from params.csv
    # library_path created for each library in the make libraries
    # function
    # extended file coming from merge reads function originally

    trimmed_file_path = os.path.join(
        library_path, '{}_trimmed.fasta'.format(run_name))

    has_adapter, has_element, has_killSequence, is_trimmedLine = [False]*4
    # set up some booleans
    trimmed_line = ''
    count, adapter_pos, adapter_dist, element_pos, element_dist = 0, -1, -1, -1, -1
    hits = []
    # booleans are declared in the outer for loop inner loop goes around for
    # each feature so if they are true for any feature booleans are
    # iteritively changed to True
    with open(trimmed_file_path, 'w') as tfp:
        for i, record in enumerate(SeqIO.parse(extended_file, 'fastq')):
            if i % 10000 == 0:
                print(i, 'records parsed')
            # keep for now convert to biopython later
            sequence = record.seq
            has_adapter, has_element, has_killSequence = False, False, False
            # not sure why reassign here?
            for feat in features:
                feature, threshold, t = feat[1].strip(
                ), feat[2].strip(), feat[3].strip()
                if t == 'remove' and feature in record.seq:
                    has_killSequence = True
                    # search for kill seq motif
                elif t == 'adapter':
                    adapter_pos, adapter_dist = kmer_search(
                        str(record.seq), feature, 6, 2)
                    adapter_len = len(feature)
                    if (adapter_pos != -1):
                        has_adapter = True
                elif t == "element":
                    element_pos, element_dist = kmer_search(
                        str(record.seq), feature, 3, 1)
                    if (element_pos != -1):
                        has_element = True
            # trim and store eligible sequences
            if not has_killSequence and has_adapter and has_element:
                adapter_end = adapter_pos+adapter_len
                trimmed_line = str(sequence)[adapter_end:element_pos]
                if len(trimmed_line) >= end_size:
                    record.letter_annotations = {}
                    record.seq = Seq(trimmed_line)  # cast back to seq object
                    hits.append(record)  # save the record
    SeqIO.write(hits, trimmed_file_path, 'fasta')
    
    # redundancy stuff I think is causing the issues look into that later

    return trimmed_file_path


def remove_duplicates(library_path, run_name, trimmed_file_path, endsize=20):
    # create rendundancy map return filtered fasta file path
    rm = RedundancyMap(trimmed_file_path, endsize)
    return rm.out


# coming soon...
