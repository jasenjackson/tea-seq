import os

from io import if_not_dir_make
from library_filter import


def make_libraries(results_dir, run_name, fastq_dir_dict, flash_exe):
    run_path = os.path.join(results_dir, run_name)
    for identifier, fastq_files in fastq_dir_dict.items:
        library_path, made_dir = if_not_dir_make(run_path, identifier)

        if made_dir:  # dir did not already exist
            combined_reads, uncombined_reads = merge_reads(
                fastq_files[0], fastq_files[1], flash_exe)

            all_reads = collate(library_path, run_name,
                                combined_reads, uncombined_reads)


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
