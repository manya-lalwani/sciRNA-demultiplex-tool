#!/usr/bin/env python3

#### imports
import sys
import gzip
from tqdm import tqdm
import difflib

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
    assert len(sequence) == len(barcode)
    count = 0
    for nuc in range(0, len(sequence)): 
        if sequence[nuc] != barcode[nuc]:
            count += 1
    return count

def one_off(sequence, barcode, barcode_length, type, diff = 1):
    '''Input: R1 sequence, expected barcode, expected length of barcode, type of barcode (rt or lig), # of nucleotides off from expected barcode to be a match (default 1)
        Output: if sequence has 1 (or diff) off match to expected barcode, returns expected barcode and start position index of match. Otherwise returns None, None
        RT: returns last match in sequence
        Lig: returns first match in sequence'''
    
    range_dict = {'rt': [len(sequence)-barcode_length, -1, -1], 'lig': [0, len(sequence)-barcode_length+1,1]}
    start, stop, step = range_dict[type]
    
    for i in range(start, stop, step):
        subseq = sequence[i:i + barcode_length]
        if how_many_off(subseq, barcode) <= diff:
            return barcode, i
            
    return None, None

def three_off(sequence, barcode, barcode_length, type, diff = 3):
    '''Input: R1 sequence, expected barcode, expected length of barcode, type of barcode (rt or lig), # of nucleotides off from expected barcode to be a match (default 1, 2, or 3)
        RT Output: If 3 off match (or less) to last two subsequences of sequence, return expected rt barcode and start position of match. Otherwise returns None, None
        Lig Output: If 3 off match (or less) to first two subsequences of sequence, return expected lig barcode and start position of match. Otherwise returns None, None'''
    
    range_map = {'rt': [len(sequence)-barcode_length, len(sequence)-barcode_length-2, -1], 'lig': [0, min(2, len(sequence)-barcode_length+1),1]}
    start, stop, step = range_map[type]

    for i in range(start, stop, step):
        subseq = sequence[i:i + barcode_length]
        if how_many_off(subseq, barcode) <= diff:
            return barcode, i
        
    return None, None

def find_closest_match(sequence, exp_barcodes, barcode_length, type, warn_threshold):
    '''Input: sequence with barcode, list of expected barcodes, barcode_length, type of barcode ("rt", "lig", or "p7"), threshold under which to add a warning to record (integer or float)
        Output: closest barcode match, boolean warn value, start index position of closest match in sequence'''
    
    closest_match, closest_match_start = None, None

    # check for 1 off differences in every set of barcode_length nucleotides, starting from the end for RT and from the beginning for ligation
    for barcode in exp_barcodes:
        closest_match, closest_match_start = one_off(sequence, barcode, barcode_length = barcode_length, type = type, diff = 1)
        if closest_match is not None:
            return closest_match, False, closest_match_start, 'one_off'
    
    # check for 2 or 3 off differences in the last 2 potential RT matches in the sequence OR
    # check for 2 or 3 off differences in the first 2 potential lig matches in the sequence
    for barcode in exp_barcodes:
        closest_match, closest_match_start = three_off(sequence, barcode, type = type, barcode_length = barcode_length, diff = 3)
        if closest_match is not None:
            return closest_match, False, closest_match_start, 'two_three_off'
        
    # remaining samples (best guess)
    closest_score = -1
    third_term = (1.3*(1)) if type == 'lig' else ((1.3*(len(sequence)-barcode_length))/len(sequence))
    max_score = (0.2*barcode_length) + (1.5*1) + third_term
    
    for barcode in exp_barcodes:
        # criteria 1: how long the match is (up to barcode_length)
        match = difflib.SequenceMatcher(None, sequence, barcode).find_longest_match(0, len(sequence), 0, len(barcode))
        match_length = match.size if match else 0

        # criteria 2: similarity of best barcode_length match, where best = least # off from barcode
        match_end = match.a + match.size
        start_val = max(0, match_end - barcode_length)
        end_val = min(len(sequence), match_end + (barcode_length - match_length)) - barcode_length + 1

        similarity = -1
        best_poss_barcode_start = 0
        for start in range(start_val, end_val):
            poss_barcode = sequence[start:start+barcode_length]
            poss_barcode_similarity = 1/(how_many_off(poss_barcode, barcode)+1)
            if poss_barcode_similarity > similarity:
                similarity = poss_barcode_similarity
                best_poss_barcode_start = start
        
        # criteria 3: distance of best poss_rt to beginning of sequence
        distance_to_start = best_poss_barcode_start/len(sequence)

        # calculate score with 3 criteria, compare to best score
        third_term = (1.3*(1/(distance_to_start+1))) if type == 'lig' else (1.3*distance_to_start)
        score = (0.2*match_length) + (1.5*similarity) + third_term

        percentage_score = (score/max_score)*100
        if percentage_score > closest_score:
            closest_score = percentage_score
            closest_match = barcode
            closest_match_start = best_poss_barcode_start
    if closest_score <= warn_threshold:
        warn_bool = True
    else:
        warn_bool = False

    return closest_match, warn_bool, closest_match_start, 'guess'

#### program functions

def check_inputs(input_fastq_folder, num_samples, num_lanes, repository_list, r1_sequence_length, umi_start_index, umi_length):
    '''Input: fastq folder with .fastq.gz files, number of samples, number of lanes, list of paths to repositories, r1 sequence length, umi start index, umi length
        Output: If fastq files are all .fastq.gz and all files are accounted for, all repositories are accounted for and plain files, and UMI does not exceed r1 sequence length, returns None. Otherwise quits program'''
    
    input_issue = False

    fastq_names = [
        f'{input_fastq_folder}/ScaleRNA_S{sample}_L{lane:03}_{file_type}_001.fastq.gz'
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
            return None  # Indicate that the file was not found

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
        sys.exit()

def compile_data(input_fastq_folder, output_fastq_folder, num_samples, num_lanes, file_type):
    '''Input: fastq folder with R1, R2, and I1 reads, fastq folder to output compiled files to, number of samples, number of lanes (up to 999), type of file (R1, R2, or I1)
        Output: returns None, writes file_type file to output fastq folder'''
    
    file_names = [
        f'{input_fastq_folder}/ScaleRNA_S{sample}_L{lane:03}_{file_type}_001.fastq.gz'
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

        if rt_start == sequence_length - rt_barcode_length:
            umis.append(sequence[umi_start:umi_start + umi_length])
            umi_qualities.append(quality[umi_start:umi_start + umi_length])
        elif rt_start > sequence_length - rt_barcode_length - umi_length:
            umis.append(sequence[umi_start:rt_start] + sequence[rt_start+rt_barcode_length:])
            umi_qualities.append(quality[umi_start:rt_start] + quality[rt_start+rt_barcode_length:])
        elif rt_start <= sequence_length - rt_barcode_length - umi_length:
            umis.append(sequence[sequence_length - umi_length:])
            umi_qualities.append(quality[sequence_length - umi_length:])

    return umis, umi_qualities

def extract_RT_barcodes(rt_barcode_repo, fastq_folder, sequence_length, warn_threshold = 55):
    '''Input: RT barcode repository file, fastq folder with R1, R2, and I1 fastq files, r1 sequence length, percentage threshold under which to warn that extracted RT barcode may be inaccurate
        Output: List of observed RT barcodes in R1 reads, list of warnings for each read (true/false), list of observed RT barcode start index for each read, RT barcode length'''

    # load in RT barcode repository, check RT barcode length of first line
    exp_rt_barcodes = [] 
    with open(rt_barcode_repo, 'r') as exp_rt_barcodes_file:
        if len((exp_rt_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_rt_barcodes_file.readline().split()[0]))
            print('Detected RT Barcode Length: ' + str(barcode_length))
        else:
            print('RT Barcode Length cannot be detected. Check RT repository for leading whitespace.')
            sys.exit()
        for line in exp_rt_barcodes_file:
            exp_rt_barcodes.append(line[0:barcode_length])

    # extract DNA sequences from R1 fastq
    r1_sequences = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'dna', sequence_length)

    # extract RT barcodes from R1 sequences and correct according to expected RT barcodes
    obs_rt_barcodes = []
    warn_rt = []
    obs_rt_start = []

    corr_pos_corr_match = 0
    wrong_pos_corr_match = 0
    one_off = 0
    two_three_off = 0
    guess = 0

    for sequence in tqdm(r1_sequences, desc = 'Extracting and correcting RT barcodes'):
        match_found = False
        if sequence[-barcode_length:] in exp_rt_barcodes:
            obs_rt_barcodes.append(sequence[-barcode_length:])
            warn_rt.append(False)
            obs_rt_start.append(len(sequence)-barcode_length)
            match_found = True
            corr_pos_corr_match += 1
        while match_found is False:
            # if barcode not in expected position, look for all barcodes throughout entire sequence
            for number, rt_barcode in enumerate(exp_rt_barcodes, 1):
                index_0 = sequence.find(rt_barcode)
                if index_0 != -1:
                    obs_rt_barcodes.append(sequence[index_0:index_0 + barcode_length])
                    warn_rt.append(False)
                    obs_rt_start.append(index_0)
                    match_found = True
                    wrong_pos_corr_match += 1
                    break
            if (number == len(exp_rt_barcodes)) and (match_found is False):
                match_found = None
        if match_found is None: # no exact match for barcode
            closest_match, warn_rt_bool, closest_match_start, method = find_closest_match(sequence, exp_rt_barcodes, type = 'rt', barcode_length = barcode_length, warn_threshold = warn_threshold) 
            obs_rt_barcodes.append(closest_match)
            warn_rt.append(warn_rt_bool)
            obs_rt_start.append(closest_match_start)
            if method == 'one_off':
                one_off += 1
            elif method == 'two_three_off':
                two_three_off += 1
            else:
                guess += 1
    
    return obs_rt_barcodes, warn_rt, obs_rt_start, barcode_length, corr_pos_corr_match, wrong_pos_corr_match, one_off, two_three_off, guess

def extract_lig_barcodes(lig_barcode_repo, fastq_folder, sequence_length, warn_threshold = 100):
    '''Input: Ligation barcode repository file, fastq folder with R1, R2, and I1 fastq files, r1 sequence length, percentage threshold under which to warn that extracted lig barcode may be inaccurate
        Output: List of observed ligation barcodes in R1 reads, list of warnings for each read (true/false), list of observed lig barcode start index for each read, lig barcode length'''
    
    # load in ligation barcode repository, check lig barcode length of first line
    exp_lig_barcodes = [] 
    with open(lig_barcode_repo, 'r') as exp_lig_barcodes_file:
        if len((exp_lig_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_lig_barcodes_file.readline().split()[0]))
            print('Detected Ligation Barcode Length: ' + str(barcode_length))
        else:
            print('RT Barcode Length cannot be detected. Check RT repository for leading whitespace.')
            sys.exit()
        for line in exp_lig_barcodes_file:
            exp_lig_barcodes.append(line[0:barcode_length])

    # extract DNA sequences from R1 fastq
    r1_sequences = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'dna', sequence_length)

    # extract ligation barcodes from R1 sequences
    obs_lig_barcodes = []
    warn_lig = []
    obs_lig_start = []

    corr_pos_corr_match = 0
    wrong_pos_corr_match = 0
    one_off = 0
    two_three_off = 0
    guess = 0

    for sequence in tqdm(r1_sequences, desc = 'Extracting ligation barcodes'):
        match_found = False
        if sequence[:barcode_length] in exp_lig_barcodes:
            obs_lig_barcodes.append(sequence[:barcode_length])
            warn_lig.append(False)
            obs_lig_start.append(0)
            corr_pos_corr_match += 1
            match_found = True  
        while match_found is False:
            # if barcode not in expected position, search for all barcodes throughout R1 sequence
            for number, lig_barcode in enumerate(exp_lig_barcodes, 1):
                index_0 = sequence.find(lig_barcode)
                if index_0 != -1:
                    obs_lig_barcodes.append(sequence[index_0:index_0 + barcode_length])
                    warn_lig.append(False)
                    obs_lig_start.append(index_0)
                    match_found = True
                    wrong_pos_corr_match += 1
                    break
            if (number == len(exp_lig_barcodes)) and (match_found is False):
                match_found = None
        if match_found is None: # no exact barcode match
            closest_match, warn_lig_bool, closest_match_start, method = find_closest_match(sequence, exp_lig_barcodes, type = 'lig', barcode_length=barcode_length, warn_threshold=warn_threshold) 
            obs_lig_barcodes.append(closest_match)
            warn_lig.append(warn_lig_bool)
            obs_lig_start.append(closest_match_start)
            if method == 'one_off':
                one_off += 1
            elif method == 'two_three_off':
                two_three_off += 1
            else:
                guess += 1
    return obs_lig_barcodes, warn_lig, obs_lig_start, barcode_length, corr_pos_corr_match, wrong_pos_corr_match, one_off, two_three_off, guess

def extract_p7_barcodes(p7_barcode_repo, fastq_folder, sequence_length):
    '''Input: P7 barcode repository file and fastq folder with R1, R2, and I1 fastq files
        Output: List of observed p7 barcodes in I1 reads, list of warnings for each read (true/false), quality scores for each observed p7 barcode'''
    
    # load in p7 barcode repository, check length of first line
    exp_p7_barcodes = [] 
    with open(p7_barcode_repo, 'r') as exp_p7_barcodes_file: # debug repeat code
        if len((exp_p7_barcodes_file.readline().split()[0])) != 0:
            barcode_length = len((exp_p7_barcodes_file.readline().split()[0]))
            print('Detected RT Barcode Length: ' + str(barcode_length))
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
    warn_p7 = []
    
    corr_match = 0
    one_two_off = 0
    guess = 0

    # if i1_sequence is 2 or less nucleotides off from expected p7 barcode, do not warn
    for sequence in tqdm(i1_sequences, desc = 'Extracting P7 barcodes'):
        if sequence[:barcode_length] in exp_p7_barcodes:
            obs_p7_barcodes.append(sequence[:barcode_length])
            warn_p7.append(False)
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
                warn_p7.append(False)
                one_two_off += 1
            else:
                warn_p7.append(True)
                guess += 1

    return obs_p7_barcodes, warn_p7, p7_qualities, corr_match, one_two_off, guess

def get_qualities(fastq_folder, obs_start, barcode_length, sequence_length):
    '''Input: fastq_folder with R1, R2, and I1 fastq files, observed start of barcode, barcode length, R1 sequence length
        Output: list of quality scores for given barcode'''

    # extract quality scores from R1 file
    r1_qualities = fastq_to_list(fastq_folder, '/R1.fastq.gz', 'quality', sequence_length=sequence_length)

    # extract barcode quality scores from R1 quality scores debug delete
    # qualities = []

    # for number, quality in enumerate(r1_qualities):
    #     qualities.append(quality[obs_start[number]:obs_start[number] + barcode_length])

    qualities = [
        quality[obs_start[number]:obs_start[number] + barcode_length]
        for number, quality in enumerate(r1_qualities)]
    
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
    combined = [rt_record + lig_record + p7_record + umi_record for rt_record, lig_record, p7_record, umi_record in zip(rt, lig, p7, umi)]
    return combined

def make_fastq(fastq_folder, type, ids, barcodes, qualities):
    '''Input: fastq folder with R1, R2, and I1 reads, type of fastq to create (Barcode, Transcripts, Warning), ids, barcodes, and qualities to input to fastq file
        Output: returns None, creates new fastq file'''

    with open(f'{fastq_folder}/{type}.fastq', 'a+') as fastq:
        for record in range(len(ids)):
            fastq.write(f'{ids[record]}{barcodes[record]}\n+\n{qualities[record]}\n')

    return None

def warn(fastq_folder, warn_rt, warn_lig, warn_p7, ids, barcodes, qualities):
    '''Input: fastq folder with R1, R2, and I1 fastq files, true/false list of warnings for each record for RT barcode, ligation, and p7, full list of ids, combined barcodes, and quality scores
        Output: returns None, creates fastq file with '''

    warn_arr = [i or j or k for i, j, k in zip(warn_rt, warn_lig, warn_p7)]
    num_unsure = sum(warn_arr)
    print('There are ' + str(num_unsure) + ' records that may be incorrect. Please check them in Warning.fastq.')

    warn_ids = []
    warn_barcodes = []
    warn_qualities = []
    warn_type = []
    for record in range(0, len(ids)): 
        if warn_arr[record]: 
            warn_ids.append(ids[record])
            warn_barcodes.append(barcodes[record])
            warn_qualities.append(qualities[record])

            warn_type_string = ''
            if warn_rt[record]:
                warn_type_string += 'RT;'
            if warn_lig[record]:
                warn_type_string += 'Lig;'
            if warn_p7[record]:
                warn_type_string += 'P7;'
            warn_type.append(warn_type_string)
        
    warn_id_and_type = [i + j for i, j in zip(warn_ids, warn_type)]
    make_fastq(fastq_folder, 'Warning', warn_id_and_type, warn_barcodes, warn_qualities)

    return None

def info(rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess, p7_corr, p7_one_two_off, p7_guess):
    '''Input: Records with correctly matched and position RT barcodes, mispositioned, one off barcodes, two or three off barcodes, and guess. Repeat for ligation barcodes. Records with correct p7 barcodes, one or two off barcodes, and guess.
        Output: None. Prints information'''
    
    print(f"The RT analysis shows {rt_corr} correct barcodes, "
          f"{rt_wrong_pos} mispositioned correct barcodes, {rt_one_off} one-off, nearby (correct place or one off) positioned barcodes, "
          f"{rt_two_three_off} two or three-off nearby barcodes, and {rt_guess} guesses (more than 3 off and possibly mispositioned).")

    print(f"The ligation analysis shows {lig_corr} correct barcodes, "
          f"{lig_wrong_pos} mispositioned correct barcodes, {lig_one_off} one-off, nearby (correct place or one off) positioned barcodes, "
          f"{lig_two_three_off} two or three-off nearby barcodes, and {lig_guess} guesses (more than 3 off and possibly mispositioned).")

    print(f"The P7 analysis shows {p7_corr} correct barcodes, "
          f"{p7_one_two_off} one or two-off barcodes, and {p7_guess} guesses (more than 2 off).")

if __name__ == "__main__":

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
    rt_warn_threshold = int(sys.argv[10])
    lig_warn_threshold = int(sys.argv[11])
    umi_start_index = int(sys.argv[12])
    umi_length = int(sys.argv[13])
    r2_sequence_length = int(sys.argv[14])
    id_character_count = int(sys.argv[15])

    #### check inputs
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
    obs_rt_barcodes, warn_rt, obs_rt_start, rt_barcode_length, rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess = extract_RT_barcodes(rt_barcode_repo, output_fastq_folder, r1_sequence_length, rt_warn_threshold)
    rt_qualities = get_qualities(output_fastq_folder, obs_rt_start, rt_barcode_length, r1_sequence_length)
    print('RT Barcodes Extracted and Corrected.')
    print('----------------------------------------------------------------------------------')
    
    #### extract observed ligation barcodes and correct them
    print('Extracting and Correcting Ligation Barcodes')
    obs_lig_barcodes, warn_lig, obs_lig_start, lig_barcode_length, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess = extract_lig_barcodes(lig_barcode_repo, output_fastq_folder, r1_sequence_length, lig_warn_threshold)
    lig_qualities = get_qualities(output_fastq_folder, obs_lig_start, lig_barcode_length, r1_sequence_length)
    print('Ligation Barcodes Extracted and Corrected.')
    print('----------------------------------------------------------------------------------')
    
    #### extract observed p7 barcodes
    print('Extracting p7 Barcodes...')
    obs_p7_barcodes, warn_p7, p7_qualities, p7_corr, p7_one_two_off, p7_guess = extract_p7_barcodes(p7_barcode_repo, output_fastq_folder, i1_sequence_length)
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

    make_fastq(output_fastq_folder, 'Barcode', ids, combined_barcodes_and_umis, combined_quality_scores)
    make_fastq(output_fastq_folder, 'Transcripts', ids, cdnas, cdna_qualities)
    warn(output_fastq_folder, warn_rt, warn_lig, warn_p7, ids, combined_barcodes_and_umis, combined_quality_scores)
    print('Result FastQ Files in ' + str(output_fastq_folder) + '.')
    print('----------------------------------------------------------------------------------')
    
    #### print information about barcode search results
    info(rt_corr, rt_wrong_pos, rt_one_off, rt_two_three_off, rt_guess, lig_corr, lig_wrong_pos, lig_one_off, lig_two_three_off, lig_guess, p7_corr, p7_one_two_off, p7_guess)
    print('----------------------------------------------------------------------------------')