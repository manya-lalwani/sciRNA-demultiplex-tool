#!/bin/bash

tool_folder="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${tool_folder}/set_params.sh"
export file_type_split=$1

split_files() {
  local input_fastq="${input_fastq_folder}/${fastq_prefix}_S1_L001_${file_type_split}_001.fastq.gz"
  # create a named pipe to avoid full decompression in memory
  mkfifo temp_pipe_${file_type_split}.fastq

  # decompress in the background using pigz (multi-threaded gzip)
  pigz -dc "${input_fastq}" > temp_pipe_${file_type_split}.fastq &  
  decompressor_pid=$!  # needed later to check if decompression is finished


  # split the decompressed stream into chunks, keeping the output compressed
  split -l $((max_records_per_chunk * 4)) temp_pipe_${file_type_split}.fastq --filter="
    chunk_dir=\"${output_fastq_folder}/Chunk_\$FILE\"
    mkdir -p \"\$chunk_dir/input_fastq_folder\" \"\$chunk_dir/output_fastq_folder\"
    output_filename=\"\${fastq_prefix}_S1_L001_${file_type_split}_001.fastq.gz\"
    pigz > \"\$chunk_dir/input_fastq_folder/\$output_filename\"
  "

  rm temp_pipe_${file_type_split}.fastq
  wait $decompressor_pid  # wait for decompression to finish
}

split_files