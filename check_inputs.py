#!/usr/bin/env python3

#### imports
import sys

#### program functions

def check_inputs(input_fastq_folder, num_samples, num_lanes, repository_list, r1_sequence_length, umi_start_index, umi_length):
    '''Input: fastq folder with .fastq.gz files, number of samples, number of lanes, list of paths to repositories, r1 sequence length, umi start index, umi length
        Output: If fastq files are all .fastq.gz and all files are accounted for, all repositories are accounted for and plain files, and UMI does not exceed r1 sequence length, returns None. Otherwise quits program'''
    
    input_issue = False

    if (num_samples != 1 or num_lanes != 1):
            print('Please only run program with chunked data on one sample and one lane at a time.')
            input_issue = True

    fastq_names = [
        f'{input_fastq_folder}/{fastq_prefix}_S{sample}_L{lane:03}_{file_type}_001.fastq.gz'
        for sample in range(1, num_samples + 1)
        for lane in range(1, num_lanes + 1)
        for file_type in ['R1', 'R2', 'I1']]

    def is_gzipped(file_name):
        try:
            with open(file_name, 'rb') as file:
                chunk = file.read(1024)
                return any(b > 127 for b in chunk)
        except FileNotFoundError:
            print(f'There is no file {file_name}')
            print('Please fix issues and try again.')
            raise ValueError('Please fix issues and try again.')

    # check that all input fastq files are present and gzipped
    for name in fastq_names:
        binary = is_gzipped(name)
        if binary is None:
            input_issue = True # had to add this
        if binary is False:
            input_issue = True
            print(f'{name} is not gunzipped. Please only input gunzipped FastQ files.')

    # check that all repositories are present and plain text
    for repository in repository_list:
        binary = is_gzipped(repository)
        if binary is True:
            input_issue = True
            print(f'{repository} is gunzipped. Please only input plain FastQ repositories.')

    # check that UMI index and length do not exceed bounds of R1 sequence length
    if umi_start_index + umi_length >= r1_sequence_length:
        input_issue = True
        print('UMI start index + UMI length is larger than R1 sequence length. Please check inputs.')
    if umi_start_index < 0:
        input_issue = True
        print('UMI start index cannot be negative. Please check inputs')

    if input_issue:
        print('Please fix issues and try again.')
        raise ValueError('Please fix issues and try again.')

if __name__ == "__main__":

    #### read in arguments
    lig_barcode_repo = sys.argv[1]
    p7_barcode_repo = sys.argv[2]
    rt_barcode_repo = sys.argv[3]

    input_fastq_folder = sys.argv[4]

    num_samples = int(sys.argv[5])
    num_lanes = int(sys.argv[6])
    r1_sequence_length = int(sys.argv[7])
    umi_start_index = int(sys.argv[8])
    umi_length = int(sys.argv[9])
    fastq_prefix = str(sys.argv[10])

    #### check inputs
    print('Checking inputs...')
    check_inputs(input_fastq_folder, num_samples, num_lanes, sys.argv[1:4], r1_sequence_length, umi_start_index, umi_length)
    print('Inputs are clean.')
    print('----------------------------------------------------------------------------------')
    
