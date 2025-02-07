#!/usr/bin/env python3

#### imports
import sys
import re
import os

#### program functions

def aggregate_metrics(output_fastq_folder, num_chunks):
    '''Aggregates metrics from each chunk file into one Aggregate_Metrics.txt file.
    Input: output fastq path, number of chunk folders
    Output: none'''

    categories = ["RT analysis", "ligation analysis", "P7 analysis"]
    summed_metrics = {category: [] for category in categories}

    text_addition = "\nRT: correct, mispositioned correct, one-off/nearby (correct or 1 place off), two or three-off/nearby, guesses (>3 off, possibly mispositioned), removed (not including removed guesses, # removed guesses = # guesses*((100-warn_threshold)/100) )\nLig: correct, mispositioned correct, one-off/nearby, two or three-off/nearby, guesses (>3 off, possibly mispositioned), removed (not including removed guesses, # guesses = same calculation as RT guesses)\nP7: correct, one or two-off, guesses (more than 2-off) (no records removed due to P7 guesses)"

    for chunk in range(1, num_chunks+1):
            chunk_metrics = f'{output_fastq_folder}/Chunk_{chunk}/output_fastq_folder/Metrics.txt'
            
            with open(chunk_metrics, 'r') as chunk_metrics_file:
                content = chunk_metrics_file.read()
                for category in categories:
                    # find lines corresponding to the category
                    match = re.search(rf"{category}.*?(?=The|$)", content, re.DOTALL)
                    if match:
                        # extract all numbers in this section
                        numbers = list(map(int, re.findall(r"(?<!\(more than )(?<!P)\b\d+\b", match.group())))
                        while len(summed_metrics[category]) < len(numbers):
                            summed_metrics[category].append(0)
                        for i, num in enumerate(numbers):
                            summed_metrics[category][i] += num

    aggregate_metrics_file = f'{output_fastq_folder}/Aggregate_Metrics.txt'
    with open(aggregate_metrics_file, 'w') as aggregate_metrics_file:
        for category, numbers in summed_metrics.items():
            aggregate_metrics_file.write(f"{category}: {', '.join(map(str, numbers))}\n")
        aggregate_metrics_file.write(text_addition)
        
    return None

def delete_extraneous_files(output_fastq_folder, num_chunks, job_name, node_count):
    '''Input: output fastq path, number of chunk files, sbatch job name, number of nodes
    Output: none
    Deletes chunk folders, sbatch out and err files, and file locks.'''

    for chunk in range(1, num_chunks+1):
        chunk_folder = f"{output_fastq_folder}/Chunk_{chunk}"
        if os.path.exists(chunk_folder):
            shutil.rmtree(chunk_folder)
    
    final_files = [f"{job_name}_split_files.out", f"{job_name}_split_files.err", f"{job_name}_final_metrics.out", f"{job_name}_final_metrics.err", f"{job_name}_change_folder_names.err", f"{job_name}_change_folder_names.out"]
    for file in final_files:
        os.remove(file)

    for node in range(1, node_count):
        chunk_processing_files = [f"{job_name}_node{node}_chunk_processing.err", f"{job_name}_node{node}_chunk_processing.out"]
        for file in chunk_processing_files:
            os.remove(file)
    
    if os.path.exists(f"{output_fastq_folder}/aggregate_data_lock.lock"):
        os.remove(f"{output_fastq_folder}/aggregate_data_lock.lock")

def convert_to_bool(bool_string):
    '''Converts different representations of true and false to python bool.'''
    if bool_string == 'false' or bool_string == 'False':
        return False
    elif bool_string == 'true' or bool_string == 'True':
        return True
    elif type(bool_string) == bool and bool_string:
        return True
    elif type(bool_string) == bool and not bool_string:
        return False
    elif bool_string == 0 or bool_string == '0':
        return False
    elif bool_string == 1 or bool_string == '1':
        return True
    else: 
        raise TypeError('Bool variable passing not working')


if __name__ == "__main__":

    #### read in arguments
    output_fastq_folder = sys.argv[1]
    num_chunks = int(sys.argv[2])
    delete_extra = convert_to_bool(sys.argv[3])
    job_name = str(sys.argv[4])
    node_count = int(sys.argv[5])

    ### aggregate metrics
    print('Aggregating Metrics...')
    aggregate_metrics(output_fastq_folder, num_chunks)
    print('Metrics aggregated')
    print('----------------------------------------------------------------------------------')

    #### delete extra files
    if delete_extra:
        print('Deleting extraneous files...')
        delete_extraneous_files(output_fastq_folder, num_chunks, job_name, node_count)
        print('Extraneous files deleted.')
        print('----------------------------------------------------------------------------------')
