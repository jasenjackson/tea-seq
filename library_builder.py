import os

from io import if_not_dir_make
from feature_search import kmer_search


def make_libraries(results_dir, run_name, fastq_dir_dict, features, flash_exe):
    run_path = os.path.join(results_dir, run_name)
    for identifier, fastq_files in fastq_dir_dict.items:
        library_path, made_dir = if_not_dir_make(run_path, identifier)
        # makes new directory in run_path

        if made_dir:  # dir did not already exist
            combined_reads, uncombined_reads = merge_reads(
                fastq_files[0], fastq_files[1], library_path, flash_exe)

            all_reads = collate(library_path, run_name,
                                combined_reads, uncombined_reads)
            trimmed_reads = feature_trim(features, library_path, run_name, combined_reads)
            



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
        cmd = '{}/./flash {} {} -d {} -M 250 --interleaved-output > {}'.format(
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
    if not os.path.exists:  # files have not been combined
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
    if os.path.exists(trimmed_file_path):
        return 0  # check if already exists
    else:
        has_adapter, has_element, has_killSequence, is_trimmedLine = [False]*4
        # set up some booleans
        trimmed_line = ''
        count, adapter_pos, adapter_dist, element_pos, element_dist = 0, -1, -1, -1, -1
        with open(trimmed_file_path) as tfp:
            with open(extended_file) as extendo:
                while True:
                    header = extendo.readline()  # read line of extendedFrags.fastq
                    sequence = extendo.readline()  # read next line fasta format prob
                    plusLine = extendo.readline()  # actuall a fastq
                    quality = extendo.readline()

                    # keep for now convert to biopython later

                    if not quality:
                        break  # end of file

                    has_adapter, has_element, has_killSequence = False, False, False
                    # not sure why reassign here?

                    for feat in features:
                        feature, threshold, t = feat[1], feat[2], feat[3]
                        if t == 'remove' and feature in sequence:
                            has_killSequence = True
                            # search for kill seq motif
                        elif t == 'adapter':
                            adapter_pos, adapter_dist = kmer_search(
                                sequence, feature, 6, 2)
                        adapter_len = len(feature)
                        if (adapter_pos != -1):
                            has_adapter = True
                        elif (type == "element"):
                            element_pos, element_dist = kmer_search(
                                sequence, feature, 3, 1)
                            if (element_pos != -1):
                                has_element = True

                    # trim and store eligible sequences
                    if ((has_killSequence == False) and (has_adapter == True) 
                        and (has_element == True)):
                        # print(count)
                        adapter_end = adapter_pos+adapter_len
                        trimmed_line = sequence[adapter_end:element_pos] + '\n'
                        if len(trimmed_line) >= end_size:
                            count += 1
                            fasta_header = ">"+header
                            new_entry = fasta_header+trimmed_line
                            trimmed_file_path.write(new_entry)
                            # need to open trimmed file path to write
    return trimmed_file_path


def remove_duplicates():
    pass
# coming soon...