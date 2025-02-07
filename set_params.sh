#!/bin/bash

#### export all variables
set -a


#### input (all fastq files) and output (empty) fastq folders
input_fastq_folder=/path/to/input_fastq_folder
output_fastq_folder=/path/to/output_fastq_folder

#### demultiplex tool
demultiplex_tool_folder=/path/to/demultiplex_tool_folder

#### barcode repository paths
lig_barcode_repo=/path/to/lig_barcode_repo
p7_barcode_repo=/path/to/p7_barcode_repo
rt_barcode_repo=/path/to/rt_barcode_repo

#### sbatch parameters (see readme)
job_name=job_name
account="account_name"
partition=partition_name
short_time_length=00:15:00
long_time_length=04:00:00
nodes=4
mail_type=ALL
mail_user=email_address@rcc.uchicago.edu
num_cpus=35 # number of CPUs per node to run on
num_cpus_file_handling=10

#### parameters (see readme)
delete_extra=False
aggregate=True
fastq_prefix='ScaleRNA'
no_scoring=0 # 0 for scoring/false, 1 for no scoring/true
chunk=true # lowercase
max_records_per_chunk=500000 # recommended 500,000
batch_size=1000 # default 1, recommend 1000 for large datasets (1-4000)

num_samples=1
num_lanes=1
r1_sequence_length=34
i1_sequence_length=10
percent_guess_to_keep=20
umi_start_index=16 # index 0, length of ligation barcode + length of linker
umi_length=8
r2_sequence_length=76
id_character_count=100 # less than 100

#### do not change
demultiplex_script=${demultiplex_tool_folder}/demultiplex_tool.py
demultiplex_script_create_chunks=${demultiplex_tool_folder}/create_chunks.sh
demultiplex_script_process_chunks=${demultiplex_tool_folder}/process_chunks.py
demultiplex_script_chunk_metrics=${demultiplex_tool_folder}/chunk_metrics.py
demultiplex_script_change_names=${demultiplex_tool_folder}/change_names.sh
demultiplex_script_check_inputs=${demultiplex_tool_folder}/check_inputs.py

NUM_NODES=1
NUM_CPUS_SHORT=1


#### stop exporting variables
set +a

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Variables set automatically."
fi