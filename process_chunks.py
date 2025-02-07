#!/usr/bin/env python3

#### imports
import sys
import subprocess
from dask import delayed, compute
from dask.distributed import Client
import ast


#### program functions
@delayed
def parallel_chunk(chunk, lig_barcode_repo, p7_barcode_repo, rt_barcode_repo, output_fastq_folder, demultiplex_script, r1_sequence_length, i1_sequence_length, percent_guess_to_keep, umi_start_index, umi_length, r2_sequence_length, id_character_count, fastq_prefix, no_scoring, batch_size, num_cpus, aggregate, delete_extra):
    '''Process one chunk using demultiplex_tool script.'''
    print(f'Processing chunk {chunk}...')

    chunk_input_fastq_folder = f'{output_fastq_folder}/Chunk_{chunk}/input_fastq_folder'
    chunk_output_fastq_folder = f'{output_fastq_folder}/Chunk_{chunk}/output_fastq_folder'
    
    arguments = [lig_barcode_repo, p7_barcode_repo, rt_barcode_repo, chunk_input_fastq_folder, chunk_output_fastq_folder, 1, 1, r1_sequence_length, i1_sequence_length, percent_guess_to_keep, umi_start_index, umi_length, r2_sequence_length, id_character_count, fastq_prefix, no_scoring, batch_size, num_cpus, True, aggregate, delete_extra]
    arguments = [str(arg) for arg in arguments]
    
    subprocess.run(['python', demultiplex_script] + arguments)

if __name__ == "__main__":

    client = Client()

    #### read in arguments
    lig_barcode_repo = sys.argv[1]
    p7_barcode_repo = sys.argv[2]
    rt_barcode_repo = sys.argv[3]

    input_fastq_folder = sys.argv[4]
    output_fastq_folder = sys.argv[5]

    num_samples = int(sys.argv[6])
    num_lanes = int(sys.argv[7])
    r1_sequence_length = int(sys.argv[8])
    i1_sequence_length = int(sys.argv[9])
    percent_guess_to_keep = int(sys.argv[10])
    umi_start_index = int(sys.argv[11])
    umi_length = int(sys.argv[12])
    r2_sequence_length = int(sys.argv[13])
    id_character_count = int(sys.argv[14])
    fastq_prefix = str(sys.argv[15])
    demultiplex_script = str(sys.argv[16])
    max_records_per_chunk = int(sys.argv[17])
    no_scoring = sys.argv[18]
    batch_size = int(sys.argv[19])
    num_cpus = int(sys.argv[20])
    node_count = int(sys.argv[21])
    chunks_per_node_list = ast.literal_eval(sys.argv[22])
    aggregate = bool(sys.argv[23])
    delete_extra = bool(sys.argv[24])

    start_chunk = 1
    for node in range(0, node_count-1):
        start_chunk += int(chunks_per_node_list[node])
    end_chunk = start_chunk + int(chunks_per_node_list[node_count-1])

    tasks = []
    for chunk in range(start_chunk, end_chunk):
        tasks.append(parallel_chunk(chunk, lig_barcode_repo, p7_barcode_repo, rt_barcode_repo, output_fastq_folder, demultiplex_script, r1_sequence_length, i1_sequence_length, percent_guess_to_keep, umi_start_index, umi_length, r2_sequence_length, id_character_count, fastq_prefix, no_scoring, batch_size, num_cpus, aggregate, delete_extra))
    compute(*tasks)

    client.close()