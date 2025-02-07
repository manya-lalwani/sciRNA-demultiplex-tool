#!/usr/bin/env python3

#### imports
import sys
import gzip
from tqdm import tqdm
import numpy as np
from optimized_functions import find_closest_match_cython
from filelock import FileLock

#### helper functions

def fastq_to_list(directory, name, type, sequence_length = 1000):
    '''Input: fastq folder, name of fastq file, type of sequence(dna, quality, or id), length of sequence to extract
        If sequence_length not given, extracts entire sequence or 1000 characters (whichever is less)
        Output: List of sequence records in file'''
    
    record_list = [] 

    sequence_type_dict = {'dna': 1, 'quality': 3, 'id': 0}
    if type in sequence_type_dict:
        begin = sequence_type_dict[type]
    else:
        print('Please specify type of extraction: dna, quality, or id')
        sys.exit()
    
    with open(directory + name, 'rb') as file: 
        file = gzip.open(file)
        for line in file.readlines()[begin::4]:
            record_list.append(line.decode('utf-8')[0:sequence_length])

    return record_list

def how_many_off(sequence, barcode):
    '''Input: sequence to compare to barcode, barcode (must be of same length)
        Output: number of nucleotides in sequence that must be changed to match barcode'''
    
    return np.sum(np.array(list(sequence)) != np.array(list(barcode)))

def process_RT_batch(batch, exp_rt_barcodes, barcode_length, no_scoring, percent_guess_to_keep):
    '''Process one dask batch of RT sequences. Check for exact matches in correct place,
    then mispositioned exact matches, then one-off matches, then two-three off matches in
    the correct place. If no_scoring is False, find best match using searching algorithm.'''

    batch_obs_rt_barcodes = []
    batch_warn_rt = []
    batch_obs_rt_start = []

    method_dict = {'corr_pos_corr_match': 0, 'wrong_pos_corr_match': 0, 'one_off': 0, 'two_three_off': 0, 'guess': 0, 'remove': 0}

    for sequence in batch:
        match_found = False
        sequence_length = len(sequence)
        sequence_end = sequence[-barcode_length:]
        if sequence_end in exp_rt_barcodes:
            batch_obs_rt_barcodes.append(sequence_end)
            batch_warn_rt.append(False)
            batch_obs_rt_start.append(sequence_length-barcode_length)
            match_found = True
            method_dict['corr_pos_corr_match'] += 1

        while match_found is False:
            # if barcode not in expected position, look for all barcodes throughout entire sequence
            for number, rt_barcode in enumerate(exp_rt_barcodes, 1):
                index_0 = sequence.find(rt_barcode)
                if index_0 != -1:
                    batch_obs_rt_barcodes.append(sequence[index_0:index_0 + barcode_length])
                    batch_warn_rt.append(False)
                    batch_obs_rt_start.append(index_0)
                    match_found = True
                    method_dict['wrong_pos_corr_match'] += 1
                    break
            if (number == len(exp_rt_barcodes)) and (match_found is False):
                match_found = None

        if not match_found: # no exact match for barcode
            closest_match, closest_score, closest_match_start, method = find_closest_match_cython(sequence, exp_rt_barcodes, barcode_type = 'rt', barcode_length = barcode_length, no_scoring = no_scoring, percent_guess_to_keep = percent_guess_to_keep) 
            batch_obs_rt_barcodes.append(closest_match)
            batch_warn_rt.append(closest_score)
            batch_obs_rt_start.append(closest_match_start)
            method_dict[method] += 1

    return batch_obs_rt_barcodes, batch_warn_rt, batch_obs_rt_start, method_dict['corr_pos_corr_match'], method_dict['wrong_pos_corr_match'], method_dict['one_off'], method_dict['two_three_off'], method_dict['guess'], method_dict['remove']

def process_lig_batch(batch, exp_lig_barcodes, barcode_length, no_scoring, percent_guess_to_keep):
    '''Process one dask batch of lig sequences. Check for exact matches in correct place,
    then mispositioned exact matches, then one-off matches, then two-three off matches in
    the correct place. If no_scoring is False, find best match using searching algorithm.'''

    batch_obs_lig_barcodes = []
    batch_warn_lig = []
    batch_obs_lig_start = []

    method_dict = {'corr_pos_corr_match': 0, 'wrong_pos_corr_match': 0, 'one_off': 0, 'two_three_off': 0, 'guess': 0, 'remove': 0}

    for sequence in batch:
        match_found = False
        sequence_length = len(sequence)
        sequence_beginning = sequence[:barcode_length]
        if sequence_beginning in exp_lig_barcodes:
            batch_obs_lig_barcodes.append(sequence_beginning)
            batch_warn_lig.append(False)
            batch_obs_lig_start.append(sequence_length-barcode_length)
            match_found = True
            method_dict['corr_pos_corr_match'] += 1
        while match_found is False:
            # if barcode not in expected position, look for all barcodes throughout entire sequence
            for number, lig_barcode in enumerate(exp_lig_barcodes, 1):
                index_0 = sequence.find(lig_barcode)
                if index_0 != -1:
                    batch_obs_lig_barcodes.append(sequence[index_0:index_0 + barcode_length])
                    batch_warn_lig.append(False)
                    batch_obs_lig_start.append(index_0)
                    match_found = True
                    method_dict['wrong_pos_corr_match'] += 1
                    break
            if (number == len(exp_lig_barcodes)) and (match_found is False):
                match_found = None
        if not match_found: # no exact match for barcode
            closest_match, closest_score, closest_match_start, method = find_closest_match_cython(sequence, exp_lig_barcodes, barcode_type = 'lig', barcode_length = barcode_length, no_scoring = no_scoring, percent_guess_to_keep = percent_guess_to_keep) 
            batch_obs_lig_barcodes.append(closest_match)
            batch_warn_lig.append(closest_score)
            batch_obs_lig_start.append(closest_match_start)
            method_dict[method] += 1

    return batch_obs_lig_barcodes, batch_warn_lig, batch_obs_lig_start, method_dict['corr_pos_corr_match'], method_dict['wrong_pos_corr_match'], method_dict['one_off'], method_dict['two_three_off'], method_dict['guess'], method_dict['remove']

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

#### program functions

def check_inputs(input_fastq_folder, num_samples, num_lanes, repository_list, r1_sequence_length, umi_start_index, umi_length):
    '''Input: fastq folder with .fastq.gz files, number of samples, number of lanes, list of paths to repositories, r1 sequence length, umi start index, umi length
        Output: If fastq files are all .fastq.gz and all files are accounted for, all repositories are accounted for and plain files, and UMI does not exceed r1 sequence length, returns None. Otherwise quits program'''
    
    input_issue = False

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
            sys.exit()

    # check that all input fastq files are present and gzipped
    for name in fastq_names:
        binary = is_gzipped(name)
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
        raise TypeError('Please fix issues and try again.')

def compile_data(input_fastq_folder, output_fastq_folder, num_samples, num_lanes, file_type):
    '''Input: fastq folder with R1, R2, and I1 reads, fastq folder to output compiled files to, number of samples, number of lanes (up to 999), type of file (R1, R2, or I1)
        Output: returns None, writes file_type file to output fastq folder'''
    
    file_names = [
        f'{input_fastq_folder}/{fastq_prefix}_S{sample}_L{lane:03}_{file_type}_001.fastq.gz'
        for sample in range(1, num_samples + 1)
        for lane in range(1, num_lanes + 1)]

    output_file_name = f'{output_fastq_folder}/{file_type}.fastq.gz'
    
    with open(output_file_name, 'ab') as output_file:
        for name in file_names:
            with open(name, 'rb') as input_file:
                output_file.write(input_file.read())

    return None

def extract_umis(fastq_folder, obs_rt_start, rt_barcode_length, sequence_length, umi_start, umi_length):
    '''Input: fastq folder with R1, R2, and I1 fastq files, RT start position, RT barcode length, R1 sequence length, UMI start position in R1 (number, index 0), UMI length
    UMI start position = length of ligation barcode + length of linker
    Output: list of UMIs in R1 reads corrected for RT barcode position, list of UMI quality scores'''

    # extract DNA sequences and quality scores from R1 fastq
    r1_sequences = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'dna', sequence_length)
    r1_qualities = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'quality', sequence_length)

    # extract UMIs from R1 sequences, correct for RT barcode mispositioning
    umis = []
    umi_qualities = []

    for number, sequence in enumerate(r1_sequences):
        rt_start = obs_rt_start[number]
        quality = r1_qualities[number]

        try:
            if rt_start == sequence_length - rt_barcode_length:
                umis.append(sequence[umi_start:umi_start + umi_length])
                umi_qualities.append(quality[umi_start:umi_start + umi_length])
            elif rt_start > sequence_length - rt_barcode_length - umi_length:
                umis.append(sequence[umi_start:rt_start] + sequence[rt_start+rt_barcode_length:])
                umi_qualities.append(quality[umi_start:rt_start] + quality[rt_start+rt_barcode_length:])
            elif rt_start <= sequence_length - rt_barcode_length - umi_length:
                umis.append(sequence[sequence_length - umi_length:])
                umi_qualities.append(quality[sequence_length - umi_length:])
        except TypeError:
            umis.append('no_scoring_to_delete')
            umi_qualities.append('no_scoring_to_delete')

    return umis, umi_qualities

def extract_p7_barcodes(p7_barcode_repo, fastq_folder, sequence_length):
    '''Input: P7 barcode repository file and fastq folder with R1, R2, and I1 fastq files
        Output: List of observed p7 barcodes in I1 reads, quality scores for each observed p7 barcode'''
    
    # load in p7 barcode repository, check length of first line
    exp_p7_barcodes = [] 
    with open(p7_barcode_repo, 'r') as exp_p7_barcodes_file:
        if len((exp_p7_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_p7_barcodes_file.readline().split()[0]))
            print('Detected P7 Barcode Length: ' + str(barcode_length))
        else:
            print('RT Barcode Length cannot be detected. Check RT repository for leading whitespace.')
            sys.exit()
        for line in exp_p7_barcodes_file:
            exp_p7_barcodes.append(line[0:barcode_length])

    # extract DNA sequences from I1 fastq
    i1_sequences = fastq_to_list(fastq_folder, '/I1.fastq.gz', 'dna', sequence_length)
    p7_qualities = fastq_to_list(fastq_folder, '/I1.fastq.gz', 'quality', sequence_length)

    # extract P7 barcodes from I1 sequences
    obs_p7_barcodes = []
    
    corr_match = 0
    one_two_off = 0
    guess = 0

    # if i1_sequence is 2 or less nucleotides off from expected p7 barcode, do not warn
    for sequence in tqdm(i1_sequences, desc = 'Extracting P7 barcodes', mininterval=5):
        if sequence[:barcode_length] in exp_p7_barcodes:
            obs_p7_barcodes.append(sequence[:barcode_length])
            corr_match += 1
        else:
            least_off = float('inf')
            obs_p7_barcode = None
            for exp_p7_barcode in exp_p7_barcodes:
                num_off = how_many_off(sequence[:barcode_length], exp_p7_barcode)
                if num_off < least_off:
                    least_off = num_off
                    obs_p7_barcode = exp_p7_barcode
            obs_p7_barcodes.append(obs_p7_barcode)
            if least_off < 3:
                one_two_off += 1
            else:
                guess += 1

    return obs_p7_barcodes, p7_qualities, corr_match, one_two_off, guess

def get_qualities(fastq_folder, obs_start, barcode_length, sequence_length):
    '''Input: fastq_folder with R1, R2, and I1 fastq files, observed start of barcode, barcode length, R1 sequence length
        Output: list of quality scores for given barcode'''

    # extract quality scores from R1 file
    r1_qualities = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'quality', sequence_length=sequence_length)

    qualities = []
    for number, quality in enumerate(r1_qualities):
        try:
            start = obs_start[number]
            end = start + barcode_length
            qualities.append(quality[start:end])
        except TypeError:
            qualities.append('no_scoring_to_delete')
    
    return qualities

def extract_cdnas(fastq_folder, sequence_length):
    '''Input: fastq folder with R1, R2, and I1 fastq files, R2 sequence length
    Output: list of cDNAs in R2 reads, list of quality scores for cDNAs'''

    # extract cDNA sequences and quality scores from R2 fastq
    cdnas = fastq_to_list(fastq_folder, '/R2.fastq.gz', 'dna', sequence_length)
    cdna_qualities = fastq_to_list(fastq_folder, '/R2.fastq.gz', 'quality', sequence_length)

    return cdnas, cdna_qualities

def combine(rt, lig, p7, umi):
    '''Input: 4 lists of strings
        Output: 1 combined list of elementwise concatenated strings'''
    combined = []
    for rt_record, lig_record, p7_record, umi_record in zip(rt, lig, p7, umi):
        if None in (rt_record, lig_record, p7_record, umi_record):
            combined.append('no_scoring_to_delete')
        else:
            combined.append(rt_record + lig_record + p7_record + umi_record)
    return combined

def make_fastq(fastq_folder, type, ids, barcodes, qualities, warn_rt, warn_lig, percent_guess_to_keep, check_for_delete = None):
    '''Input: fastq folder with R1, R2, and I1 reads, type of fastq to create (Barcode, Transcripts, Warning), ids, barcodes, and qualities to input to fastq file
        Output: returns None, creates new fastq file'''

    warn_rt_scores = [x for x in warn_rt if x is not False]
    warn_lig_scores = [x for x in warn_lig if x is not False]
    if warn_rt_scores:
        rt_threshold = np.percentile(warn_rt_scores, 100 - percent_guess_to_keep)
    else:
        rt_threshold = -1
    if warn_lig_scores:
        lig_threshold = np.percentile(warn_lig_scores, 100 - percent_guess_to_keep)
    else:
        lig_threshold = -1
    
    
    with open(f'{fastq_folder}/{type}.fastq', 'a+') as fastq:
        if check_for_delete:
            for record in range(len(ids)):
                if 'no_scoring_to_delete' not in check_for_delete[record]:
                    if (warn_rt[record] is False or warn_rt[record] > rt_threshold) and (warn_lig[record] is False or warn_lig[record] > lig_threshold):
                        fastq.write(f"{ids[record]}{barcodes[record].strip()}\n+\n{qualities[record].strip()}\n")
        else:
            for record in range(len(ids)):
                if 'no_scoring_to_delete' not in qualities[record]:
                    if (warn_rt[record] is False or warn_rt[record] > rt_threshold) and (warn_lig[record] is False or warn_lig[record] > lig_threshold):
                        fastq.write(f"{ids[record]}{barcodes[record].strip()}\n+\n{qualities[record].strip()}\n")


    return None

def make_aggregate_fastq(fastq_folder, ids, combined_barcodes_and_umis, combined_quality_scores, cdnas, cdna_qualities, warn_rt, warn_lig, percent_guess_to_keep):
    '''Combines chunk fastq files into two aggregate Barcode and Transcripts fastq. Returns none.'''

    outer_output_folder = fastq_folder.rsplit('/', 2)[0]
    barcode_file = f'{outer_output_folder}/Barcode.fastq.gz'
    transcript_file = f'{outer_output_folder}/Transcripts.fastq.gz'

    barcode_data = []
    transcript_data = []

    warn_rt_scores = [x for x in warn_rt if x is not False]
    warn_lig_scores = [x for x in warn_lig if x is not False]
    if warn_rt_scores:
        rt_threshold = np.percentile(warn_rt_scores, 100 - percent_guess_to_keep)
    else:
        rt_threshold = -1
    if warn_lig_scores:
        lig_threshold = np.percentile(warn_lig_scores, 100 - percent_guess_to_keep)
    else:
        lig_threshold = -1

    for record in range(len(ids)):
        if 'no_scoring_to_delete' not in combined_quality_scores[record]:
            if (warn_rt[record] is False or warn_rt[record] > rt_threshold) and (warn_lig[record] is False or warn_lig[record] > lig_threshold):
                barcode_data.append(f'{ids[record]}{combined_barcodes_and_umis[record]}\n+\n{combined_quality_scores[record]}\n')
                transcript_data.append(f"{ids[record]}{cdnas[record].strip()}\n+\n{cdna_qualities[record].strip()}\n")

    barcode_bytes = ''.join(barcode_data).encode('utf-8')
    transcript_bytes = ''.join(transcript_data).encode('utf-8')

    aggregate_data_lock = FileLock(f"{outer_output_folder}/aggregate_data_lock.lock")

    with aggregate_data_lock:
        with gzip.open(barcode_file, 'ab') as barcode_output, gzip.open(transcript_file, 'ab') as transcript_output:
            barcode_output.write(barcode_bytes)
            transcript_output.write(transcript_bytes)

    return None

def metrics(rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess, rt_remove, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess, lig_remove, p7_corr, p7_one_two_off, p7_guess, output_fastq_folder, percent_guess_to_keep):
    '''Input: Records with correctly matched and position RT barcodes, mispositioned, one off barcodes, two or three off barcodes, and guess. Repeat for ligation barcodes. Records with correct p7 barcodes, one or two off barcodes, and guess.
        Output: None. Prints information'''
    
    metrics = (f"The RT analysis shows {rt_corr} correct barcodes, "
            f"{rt_wrong_pos} mispositioned correct barcodes, {rt_one_off} one-off, nearby (correct place or one off) positioned barcodes, "
            f"{rt_two_three_off} two or three-off nearby barcodes, and {rt_guess} guesses (more than 3 off and possibly mispositioned). {rt_remove} records "
            f"have been removed due to failure to find RT barcode with basic methods (if nonzero, no_scoring = True) \n\n"

            f"The ligation analysis shows {lig_corr} correct barcodes, "
            f"{lig_wrong_pos} mispositioned correct barcodes, {lig_one_off} one-off, nearby (correct place or one off) positioned barcodes, "
            f"{lig_two_three_off} two or three-off nearby barcodes, and {lig_guess} guesses (more than 3 off and possibly mispositioned). {lig_remove} records "
            f"have been removed due to failure to find ligation barcode with basic methods (if nonzero, no_scoring = True) \n\n"

            f"The P7 analysis shows {p7_corr} correct barcodes, "
            f"{p7_one_two_off} one or two-off barcodes, and {p7_guess} guesses (more than 2 off).\n")
    
    print(metrics)

    with open(f'{output_fastq_folder}/Metrics.txt', 'w') as metrics_file:
        metrics_file.write(metrics)

def extract_RT_barcodes_parallel(rt_barcode_repo, fastq_folder, sequence_length, no_scoring, batch_size = 1, percent_guess_to_keep = 0, num_cpus = 1):
    '''Input: RT barcode repository file, fastq folder with R1, R2, and I1 fastq files, r1 sequence length, percentage threshold under which to warn that extracted RT barcode may be inaccurate
        Output: List of observed RT barcodes in R1 reads, list of scores for each read, list of observed RT barcode start index for each read, RT barcode length'''

    # load in RT barcode repository, check RT barcode length of first line
    exp_rt_barcodes = set()
    with open(rt_barcode_repo, 'r') as exp_rt_barcodes_file:
        if len((exp_rt_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_rt_barcodes_file.readline().split()[0]))
            print('Detected RT Barcode Length: ' + str(barcode_length))
        else:
            print('RT Barcode Length cannot be detected. Check RT repository for leading whitespace.')
            sys.exit()
        for line in exp_rt_barcodes_file:
            exp_rt_barcodes.add(line[0:barcode_length])

    # extract DNA sequences from R1 fastq
    print('Extracting R1 sequences...')
    r1_sequences = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'dna', sequence_length)

    # Split sequences into batches
    print('Splitting into batches...')
    batches = [r1_sequences[i:i + batch_size] for i in range(0, len(r1_sequences), batch_size)]

    obs_rt_barcodes = []
    warn_rt = []
    obs_rt_start = []
    corr_pos_corr_match = 0
    wrong_pos_corr_match = 0
    one_off = 0
    two_three_off = 0
    guess = 0
    remove = 0
    
    for batch in tqdm(batches, "Extracting RT Barcodes...", mininterval=10):
        batch_obs_rt_barcodes, batch_warn_rt, batch_obs_rt_start, batch_corr_pos_corr_match, batch_wrong_pos_corr_match, batch_one_off, batch_two_three_off, batch_guess, batch_remove = process_RT_batch(batch, exp_rt_barcodes, barcode_length, no_scoring, percent_guess_to_keep)
        obs_rt_barcodes.extend(batch_obs_rt_barcodes)
        warn_rt.extend(batch_warn_rt)
        obs_rt_start.extend(batch_obs_rt_start)
        corr_pos_corr_match += batch_corr_pos_corr_match
        wrong_pos_corr_match += batch_wrong_pos_corr_match
        one_off += batch_one_off
        two_three_off += batch_two_three_off
        guess += batch_guess
        remove += batch_remove

    return obs_rt_barcodes, warn_rt, obs_rt_start, barcode_length, corr_pos_corr_match, wrong_pos_corr_match, one_off, two_three_off, guess, remove, batches

def extract_lig_barcodes_parallel(batches, lig_barcode_repo, no_scoring, fastq_folder, sequence_length, percent_guess_to_keep = 0, num_cpus = 1):
    '''Input: Ligation barcode repository file, fastq folder with R1, R2, and I1 fastq files, r1 sequence length, percentage threshold under which to warn that extracted lig barcode may be inaccurate
        Output: List of observed ligation barcodes in R1 reads, list of scores for each read, list of observed lig barcode start index for each read, lig barcode length'''
    
    # load in ligation barcode repository, check lig barcode length of first line
    exp_lig_barcodes = set()
    with open(lig_barcode_repo, 'r') as exp_lig_barcodes_file:
        if len((exp_lig_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_lig_barcodes_file.readline().split()[0]))
            print('Detected Ligation Barcode Length: ' + str(barcode_length))
        else:
            print('Ligation Barcode Length cannot be detected. Check RT repository for leading whitespace.')
            sys.exit()
        for line in exp_lig_barcodes_file:
            exp_lig_barcodes.add(line[0:barcode_length])
    
    obs_lig_barcodes = []
    warn_lig = []
    obs_lig_start = []

    corr_pos_corr_match = 0
    wrong_pos_corr_match = 0
    one_off = 0
    two_three_off = 0
    guess = 0
    remove = 0

    for batch in tqdm(batches, "Extracting Lig Barcodes...", mininterval=10):
        batch_obs_lig_barcodes, batch_warn_lig, batch_obs_lig_start, batch_corr_pos_corr_match, batch_wrong_pos_corr_match, batch_one_off, batch_two_three_off, batch_guess, batch_remove = process_lig_batch(batch, exp_lig_barcodes, barcode_length, no_scoring, percent_guess_to_keep)
        obs_lig_barcodes.extend(batch_obs_lig_barcodes)
        warn_lig.extend(batch_warn_lig)
        obs_lig_start.extend(batch_obs_lig_start)
        corr_pos_corr_match += batch_corr_pos_corr_match
        wrong_pos_corr_match += batch_wrong_pos_corr_match
        one_off += batch_one_off
        two_three_off += batch_two_three_off
        guess += batch_guess
        remove += batch_remove

    return obs_lig_barcodes, warn_lig, obs_lig_start, barcode_length, corr_pos_corr_match, wrong_pos_corr_match, one_off, two_three_off, guess, remove

if __name__ == "__main__":

    #### read in arguments
    lig_barcode_repo = str(sys.argv[1])
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
    no_scoring = convert_to_bool(sys.argv[16])
    batch_size = int(sys.argv[17])
    num_cpus = int(sys.argv[18])
    chunk = convert_to_bool(sys.argv[19])
    aggregate = convert_to_bool(sys.argv[20])
    delete_extra = convert_to_bool(sys.argv[21])

    ## check inputs
    print('Checking inputs...')
    check_inputs(input_fastq_folder, num_samples, num_lanes, sys.argv[1:4], r1_sequence_length, umi_start_index, umi_length)
    print('Inputs are clean.')
    print('----------------------------------------------------------------------------------')

    #### compile data to get all R1/R2/I1 data together
    print('Compiling FastQs into 3 Files (R1, R2, I2)...')
    compile_data(input_fastq_folder, output_fastq_folder, num_samples, num_lanes, file_type = 'R1')
    compile_data(input_fastq_folder, output_fastq_folder, num_samples, num_lanes, file_type = 'R2')
    compile_data(input_fastq_folder, output_fastq_folder, num_samples, num_lanes, file_type = 'I1')
    print('FastQs combined. New files in ' + str(output_fastq_folder) + '.')
    print('----------------------------------------------------------------------------------')

    #### extract observed RT barcodes and correct them according to repository
    print('Extracting and Correcting RT Barcodes...')
    obs_rt_barcodes, warn_rt, obs_rt_start, rt_barcode_length, rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess, rt_remove, r1_batches = extract_RT_barcodes_parallel(rt_barcode_repo, output_fastq_folder, r1_sequence_length, no_scoring, batch_size, percent_guess_to_keep, num_cpus)
    rt_qualities = get_qualities(output_fastq_folder, obs_rt_start, rt_barcode_length, r1_sequence_length)
    print('RT Barcodes Extracted and Corrected.')
    print('----------------------------------------------------------------------------------')


    #### extract observed ligation barcodes and correct them
    print('Extracting and Correcting Ligation Barcodes')
    obs_lig_barcodes, warn_lig, obs_lig_start, lig_barcode_length, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess, lig_remove = extract_lig_barcodes_parallel(r1_batches, lig_barcode_repo, no_scoring, output_fastq_folder, r1_sequence_length, percent_guess_to_keep, num_cpus)
    lig_qualities = get_qualities(output_fastq_folder, obs_lig_start, lig_barcode_length, r1_sequence_length)
    print('Ligation Barcodes Extracted and Corrected.')
    print('----------------------------------------------------------------------------------')

    
    #### extract observed p7 barcodes
    print('Extracting p7 Barcodes...')
    obs_p7_barcodes, p7_qualities, p7_corr, p7_one_two_off, p7_guess = extract_p7_barcodes(p7_barcode_repo, output_fastq_folder, i1_sequence_length)
    print('P7 Barcodes Extracted.')
    print('----------------------------------------------------------------------------------')

    #### extract UMIs (correct for mispositioned RT barcode)
    print('Extracting and correcting UMIs...')
    umis, umi_qualities = extract_umis(output_fastq_folder, obs_rt_start, rt_barcode_length, r1_sequence_length, umi_start_index, umi_length)
    print('UMIs Extracted and Corrected.')
    print('----------------------------------------------------------------------------------')
    
    #### extract cDNA
    print('Extracting cDNAs...')
    cdnas, cdna_qualities = extract_cdnas(output_fastq_folder, r2_sequence_length)
    print('CDNAs extracted.')
    print('----------------------------------------------------------------------------------')

    #### create new fastq files
    print('Creating new FastQ Files for barcodes, transcripts, and warnings (possible incorrect processing)...')
    ids = fastq_to_list(output_fastq_folder, '/R1.fastq.gz', 'id', id_character_count)
    combined_barcodes_and_umis = combine(obs_rt_barcodes, obs_lig_barcodes, obs_p7_barcodes, umis)
    combined_quality_scores = combine(rt_qualities, lig_qualities, p7_qualities, umi_qualities)

    make_fastq(output_fastq_folder, 'Barcode', ids, combined_barcodes_and_umis, combined_quality_scores, warn_rt, warn_lig, percent_guess_to_keep)
    make_fastq(output_fastq_folder, 'Transcripts', ids, cdnas, cdna_qualities, warn_rt, warn_lig, percent_guess_to_keep, check_for_delete = combined_quality_scores)
    print('Result FastQ Files in ' + str(output_fastq_folder) + '.')
    print('----------------------------------------------------------------------------------')

    #### create aggregated fastq files
    if (chunk and (aggregate or delete_extra)):
        make_aggregate_fastq(output_fastq_folder, ids, combined_barcodes_and_umis, combined_quality_scores, cdnas, cdna_qualities, warn_rt, warn_lig, percent_guess_to_keep)
        print('Aggregated FastQ Files in ' + str(output_fastq_folder.rsplit('/', 2)[0]) + '.')
        print('----------------------------------------------------------------------------------')

    #### print metrics of barcode search results and saves to Metrics.txt
    metrics(rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess, rt_remove, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess, lig_remove, p7_corr, p7_one_two_off, p7_guess, output_fastq_folder, percent_guess_to_keep)
    print('----------------------------------------------------------------------------------')

