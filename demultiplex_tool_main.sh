#!/bin/bash

#### import variables
tool_folder="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${tool_folder}/set_params.sh"


#### helper functions

calculate_chunks_per_node() {
    local num_chunks=$1
    local num_nodes=$2
    local chunks_per_node=$((num_chunks / num_nodes))
    local remainder=$((num_chunks % num_nodes))
    
    local chunks_per_node_list=()
    for ((i=0; i<num_nodes; i++)); do
        if ((i < remainder)); then
            chunks_per_node_list+=($((chunks_per_node + 1)))
        else
            chunks_per_node_list+=($chunks_per_node)
        fi
    done
    printf "[%s]\n" "$(IFS=', '; echo "${chunks_per_node_list[*]}")"
}

submit_sbatch() {
    local job_dependency=$1
    local job_name=$2
    local output=${2}.out
    local error=${2}.err
    local partition=$3
    local account=$4
    local nodes=$5
    local ntasks_per_node=$6
    local time=$7
    local mail_type=${8}
    local mail_user=${9}
    local script_type=${10}
    shift 10

    local command="$1"
    shift
    local args="$@"

    if [ $script_type == 'bash' ]; then
        if [ -z "$job_dependency" ]; then
            job_id=$(sbatch --parsable <<EOT
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$output
#SBATCH --error=$error
#SBATCH --partition=$partition
#SBATCH --account=$account
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$ntasks_per_node
#SBATCH --time=$time
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$mail_user
bash $command $args
EOT
        )
        else
            job_id=$(sbatch --parsable <<EOT
#!/bin/bash
#SBATCH --dependency=$job_dependency
#SBATCH --job-name=$job_name
#SBATCH --output=$output
#SBATCH --error=$error
#SBATCH --partition=$partition
#SBATCH --account=$account
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$ntasks_per_node
#SBATCH --time=$time
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$mail_user
bash $command $args
EOT
        )
        fi
    else
        if [ -z "$job_dependency" ]; then
            job_id=$(sbatch --parsable <<EOT
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$output
#SBATCH --error=$error
#SBATCH --partition=$partition
#SBATCH --account=$account
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$ntasks_per_node
#SBATCH --time=$time
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$mail_user
python3 $command $args
EOT
        )
        else
            job_id=$(sbatch --parsable <<EOT
#!/bin/bash
#SBATCH --dependency=$job_dependency
#SBATCH --job-name=$job_name
#SBATCH --output=$output
#SBATCH --error=$error
#SBATCH --partition=$partition
#SBATCH --account=$account
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$ntasks_per_node
#SBATCH --time=$time
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$mail_user
python3 $command $args
EOT
        )
        fi
    fi
    echo $job_id
}





#### run program
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

    cd $output_fastq_folder
    
    #### if short file, submit a job to run demultiplex_tool.py
    if [ "$chunk" = false ]; then
        submit_sbatch "" $job_name $partition $account $NUM_NODES $NUM_CPUS_SHORT $short_time_length $mail_type $mail_user "python" $demultiplex_script $lig_barcode_repo $p7_barcode_repo $rt_barcode_repo $input_fastq_folder $output_fastq_folder $num_samples $num_lanes $r1_sequence_length $i1_sequence_length $percent_guess_to_keep $umi_start_index $umi_length $r2_sequence_length $id_character_count $fastq_prefix $no_scoring $batch_size $num_cpus $chunk $aggregate $delete_extra
    
    #### if long file, submit a job to split the file up into chunks
    else
        #### check inputs
        python $demultiplex_script_check_inputs $lig_barcode_repo $p7_barcode_repo $rt_barcode_repo $input_fastq_folder $num_samples $num_lanes $r1_sequence_length $umi_start_index $umi_length $fastq_prefix

        #### split file
        file_types=('R1' 'R2' 'I1')
        all_job_ids=""
        for file_type in "${file_types[@]}"; do
            job_id=$(submit_sbatch "" "${job_name}_split_files" $partition $account $NUM_NODES $num_cpus_file_handling $long_time_length $mail_type $mail_user "bash" $demultiplex_script_create_chunks $file_type)
            if [ -z "$all_job_ids" ]; then
                all_job_ids="afterok:$job_id"
            else
                all_job_ids="$all_job_ids:$job_id"
            fi
        done

        num_chunks=-1
        submit_sbatch "$all_job_ids" "${job_name}_change_folder_names" $partition $account $NUM_NODES $NUM_CPUS_SHORT $short_time_length $mail_type $mail_user "bash" $demultiplex_script_change_names

        #### get out number of chunks and calculate chunks per node
        sleep 1000 # for long files
        temp_file=false
        while [[ "$temp_file" == false ]]; do
            if [[ -f "${output_fastq_folder}/num_chunks.txt" ]]; then
                num_chunks=$(cat "${output_fastq_folder}/num_chunks.txt")
                chunks_per_node_list=$(calculate_chunks_per_node $num_chunks $nodes)
                rm -f "${output_fastq_folder}/num_chunks.txt"
                temp_file=true
            else
                sleep 60
            fi
        done
        echo FastQs split. New files in $output_fastq_folder.
        echo There are $num_chunks chunks to be processed.
        echo ----------------------------------------------------------------------------------

        #### submit individual jobs for every node to process chunks
        chunks_per_node_stripped=$(echo "$chunks_per_node_list" | tr -d '[]')
        IFS=', ' read -r -a chunks_per_node_array <<< "$chunks_per_node_stripped"
        
        node_count=1
        all_job_ids=""
        for alloc_chunks in "${chunks_per_node_array[@]}"; do
            job_id=$(submit_sbatch "" "${job_name}_node${node_count}_chunk_processing" $partition $account $NUM_NODES $num_cpus $long_time_length $mail_type $mail_user "python" $demultiplex_script_process_chunks $lig_barcode_repo $p7_barcode_repo $rt_barcode_repo $input_fastq_folder $output_fastq_folder $num_samples $num_lanes $r1_sequence_length $i1_sequence_length $percent_guess_to_keep $umi_start_index $umi_length $r2_sequence_length $id_character_count $fastq_prefix $demultiplex_script $max_records_per_chunk $no_scoring $batch_size $num_cpus $node_count $chunks_per_node_list $aggregate $delete_extra)
            
            if [ -z "$all_job_ids" ]; then
                all_job_ids="afterok:$job_id"
            else
                all_job_ids="$all_job_ids:$job_id"
            fi
            ((node_count++))
        done

        #### aggregate metrics from all chunks
        submit_sbatch "$all_job_ids" "${job_name}_final_metrics" $partition $account $NUM_NODES $NUM_CPUS_SHORT $short_time_length $mail_type $mail_user "python" \
            $demultiplex_script_chunk_metrics \
            $output_fastq_folder \
            $num_chunks \
            $delete_extra \
            $job_name \
            $node_count 
    fi
    echo All processes running.
fi
