#### barcode repository paths
lig_barcode_repo=/project/zhuzhuzhang/mlalwani/fastqs/test_files/lig_barcodes.gz
p7_barcode_repo=/project/zhuzhuzhang/mlalwani/fastqs/test_files/p7_barcodes.gz
rt_barcode_repo=/project/zhuzhuzhang/mlalwani/fastqs/test_files/rt_barcodes.txt

#### input (all fastq files) and output (empty) fastq folders
input_fastq_folder=/project/zhuzhuzhang/mlalwani/fastqs/input_fastq_folder
output_fastq_folder=/project/zhuzhuzhang/mlalwani/fastqs/test_files/output_fastqs

#### parameters
num_samples=1
num_lanes=4
r1_sequence_length=34
i1_sequence_length=10
rt_warn_threshold=55
lig_warn_threshold=100
umi_start_index=16 # index 0, length of ligation barcode + length of linker
umi_length=8
r2_sequence_length=76
id_character_count=100 # less than 100


#### run demultiplex tool
demultiplex_script=/project/zhuzhuzhang/mlalwani/manya_demultiplex_tool/demultiplex_tool.py

python3 $demultiplex_script \
    $lig_barcode_repo \
    $p7_barcode_repo \
    $rt_barcode_repo \
    $input_fastq_folder \
    $output_fastq_folder \
    $num_samples \
    $num_lanes \
    $r1_sequence_length \
    $i1_sequence_length \
    $rt_warn_threshold \
    $lig_warn_threshold \
    $umi_start_index \
    $umi_length \
    $r2_sequence_length \
    $id_character_count 
