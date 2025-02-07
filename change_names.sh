#!/bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

    chunk_number_file=1
    for chunk_dir in Chunk_*; do
        if [[ -d "$chunk_dir" ]]; then
            new_chunk_name="Chunk_${chunk_number_file}"

            if [[ "$chunk_dir" != "$new_chunk_name" ]]; then
                mv "$chunk_dir" "$new_chunk_name"
            fi
            ((chunk_number_file++))
        fi
    done
    echo "$((chunk_number_file - 1))" > "${output_fastq_folder}/num_chunks.txt"
fi

