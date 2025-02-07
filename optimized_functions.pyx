# optimized_functions.pyx

import difflib
from libc.string cimport memcpy
from libc.stdlib cimport malloc, free
import numpy as np
from cython.parallel cimport parallel, prange

# helper function for counting mismatches
cpdef int count_mismatches(char* seq_ptr, char* barcode_ptr, int barcode_length, int start_index):
    cdef int off_count = 0
    for j in range(barcode_length):
        if seq_ptr[start_index + j] != barcode_ptr[j]:
            off_count += 1
            if off_count > 3:  # don't check further if mismatches exceed 3
                return off_count
    return off_count

# calculate how many positions are different between two strings
cpdef int how_many_off(str s1, str s2):
    cdef int count = 0
    cdef int i, len_s1
    len_s1 = len(s1)
    for i in range(len_s1):
        if s1[i] != s2[i]:
            count += 1
    return count

cpdef tuple find_closest_match_cython(
    str sequence, 
    set exp_barcodes, 
    int barcode_length, 
    str barcode_type, 
    bint no_scoring,  
    float percent_guess_to_keep
):
    cdef int i, seq_len, start, stop, step, off_count, match_length, match_end, start_val, end_val
    cdef int best_poss_barcode_start
    cdef float score, third_term, percentage_score, closest_score
    cdef str closest_match, mismatch_type
    cdef char* seq_ptr
    cdef char* bar_ptr
    cdef double similarity
    
    # convert sequence to C string for efficiency
    seq_len = len(sequence)
    seq_ptr = <char*>malloc(seq_len + 1)  # +1 for null terminator
    if seq_ptr == NULL:
        raise MemoryError("Cannot allocate memory for seq_ptr")
    
    cdef bytes encoded_sequence = sequence.encode('ascii')
    for i in range(seq_len):
        seq_ptr[i] = encoded_sequence[i]
    seq_ptr[seq_len] = 0  # null-terminate the C string

    closest_match = None
    closest_match_start = -1
    mismatch_type = ""
    closest_score = -1

    if barcode_type == 'rt':
        start, stop, step = seq_len - barcode_length, -1, -1
    elif barcode_type == 'lig':
        start, stop, step = 0, seq_len - barcode_length + 1, 1

    # Check for 1 off differences
    cdef bytes encoded_barcode 
    for barcode in exp_barcodes:
        bar_ptr = <char*>malloc(barcode_length + 1)  # allocate memory for barcode
        if bar_ptr == NULL:
            free(seq_ptr)
            raise MemoryError("Failed to allocate memory for bar_ptr")
        
        encoded_barcode = barcode.encode('ascii')
        for i in range(barcode_length):
            bar_ptr[i] = encoded_barcode[i]
        bar_ptr[barcode_length] = 0  # null-terminate the barcode pointer

        for i in range(start, stop, step):
            off_count = count_mismatches(seq_ptr, bar_ptr, barcode_length, i)
            if off_count <= 1:
                closest_match = barcode
                closest_match_start = i
                mismatch_type = 'one_off'
                free(bar_ptr)  # free the memory before returning
                free(seq_ptr)
                return closest_match, False, closest_match_start, mismatch_type

        free(bar_ptr)  # free memory before moving to the next barcode

    # Check for 2 or 3 off differences in specific locations based on barcode type
    cdef int two_three_start, two_three_stop, two_three_step
    if barcode_type == 'rt':
        two_three_start, two_three_stop, two_three_step = seq_len - barcode_length, seq_len - barcode_length - 2, -1
    elif barcode_type == 'lig':
        two_three_start, two_three_stop, two_three_step = 0, min(2, seq_len - barcode_length + 1), 1

    for barcode in exp_barcodes:
        bar_ptr = <char*>malloc(barcode_length + 1)

        if bar_ptr == NULL:
            free(seq_ptr)
            raise MemoryError("Failed to allocate memory for bar_ptr")

        encoded_barcode = barcode.encode('ascii')
        for i in range(barcode_length):
            bar_ptr[i] = encoded_barcode[i]
        bar_ptr[barcode_length] = 0  # null-terminate the barcode pointer

        for i in range(two_three_start, two_three_stop, two_three_step):
            off_count = count_mismatches(seq_ptr, bar_ptr, barcode_length, i)
            if off_count <= 3:
                closest_match = barcode
                closest_match_start = i
                mismatch_type = 'two_three_off'
                free(bar_ptr)
                free(seq_ptr)
                return closest_match, False, closest_match_start, mismatch_type

        free(bar_ptr)  # free memory before moving to the next barcode

    # Free seq_ptr memory at the end of the function
    free(seq_ptr)

    if no_scoring is True:
        return None, False, None, 'removed'

    # Scoring based on match length, similarity, and distance
    max_score = (0.2 * barcode_length) + (1.5 * 1) + (1.3 * (1) if barcode_type == 'lig' else (1.3 * (len(sequence) - barcode_length) / len(sequence)))

    for barcode in exp_barcodes:
        # Criteria 1: Find longest match
        match = difflib.SequenceMatcher(None, sequence, barcode).find_longest_match(0, len(sequence), 0, len(barcode))
        match_length = match.size if match else 0

        # Criteria 2: Best match similarity
        match_end = match.a + match.size
        start_val = max(0, match_end - barcode_length)
        end_val = min(len(sequence), match_end + (barcode_length - match_length)) - barcode_length + 1

        similarity = -1
        best_poss_barcode_start = 0
        for start in range(start_val, end_val):
            poss_barcode = sequence[start:start+barcode_length]
            poss_barcode_similarity = 1 / (how_many_off(poss_barcode, barcode) + 1)
            if poss_barcode_similarity > similarity:
                similarity = poss_barcode_similarity
                best_poss_barcode_start = start

        # Criteria 3: Distance of best match to the start
        distance_to_start = best_poss_barcode_start / len(sequence)

        third_term = (1.3 * (1 / (distance_to_start + 1))) if barcode_type == 'lig' else (1.3 * distance_to_start)
        score = (0.2 * match_length) + (1.5 * similarity) + third_term

        # Calculate percentage score
        percentage_score = (score / max_score) * 100
        if percentage_score > closest_score:
            closest_score = percentage_score
            closest_match = barcode
            closest_match_start = best_poss_barcode_start

    return closest_match, closest_score, closest_match_start, 'guess'