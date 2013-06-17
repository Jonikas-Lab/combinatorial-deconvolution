#! /usr/bin/env python

"""
Basic functions for combinatorial deconvolution of mutant pools, i.e. matching observed deepseq readcounts for each insertion in each pool to the expected presence/absence codewords of each original sample in each pool, to determine which sequenced insertion corresponds to which original sample.

Basic experimental pipeline (and nomenclature):
 1) Biological SAMPLES are pooled into POOLS, based on predetermined CODEWORDS, which are then used as "EXPECTED" CODEWORDS.
    (sometimes the "samples" are already pools themselves, and the "pools" are called superpools, but that doesn't matter here)
    This is all done in the separate combinatorial_pooling package.
 2) After each POOL is sequenced, each INSERTION found in any pool will have an "OBSERVED" CODEWORD based on 
    which pools it was present in.  Those "observed" codewords are then matched to the closest "expected" codeword, 
    to determine which INSERTION corresponds to which original SAMPLE - this is what this module is for. 

 -- Weronika Patena, 2013
"""

# standard library
from __future__ import division
import os
import collections
import unittest
import random
import math
from collections import defaultdict
from string import ascii_uppercase # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', for well numbers
# other packages
import numpy
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import general_utilities, plotting_utilities
import binary_code_utilities
import mutant_analysis_classes

class DeconvolutionError(Exception):
    """ Exceptions in the deconvolution_utilities module."""
    pass


################################################## Various utilities #######################################################

def get_original_384_well_numbers(sample_number, zero_padding=True):
    """ Get a well name (A1-P24) form a sample number (1-384) assuming the 384-well plate was split into 4 staggered 96-well plates.
    
    It's a complicated setup!  The assumption is that the original 384-well plate was split into four staggered 96-well plates, 
     with every other row/column being in the same plate, like this:
                col1:       col2:       col3:       col4:       col5:       ...     col23:      col24:
        rowA:   plate1-A1   plate2-A1   plate1-A2   plate2-A2   plate1-A3   ...     plate1-A12  plate2-A12
        rowB:   plate3-A1   plate4-A1   plate3-A2   plate4-A2   plate3-A3   ...     plate3-A12  plate4-A12
        rowC:   plate1-B1   plate2-B1   plate1-B2   plate2-B2   plate1-B3   ...     plate1-B12  plate2-B12
                ...         ...         ...         ...         ...         ...     ...         ...
        rowO:   plate1-H1   plate2-H1   plate1-H2   plate2-H2   plate1-H3   ...     plate1-H12  plate2-H12
        rowP:   plate3-H1   plate4-H1   plate3-H2   plate4-H2   plate3-H3   ...     plate3-H12  plate4-H12

     And then that the resulting four 96-well plates were taken in order, with each well in order, 
      so the final sample order would be like this: plate1-A1, plate1-A2, ..., plate2-A1, ..., plate4-H12.

     So the full conversion from final number to original 384-well number is like this:
         - sample 1 would be plate1-A1, i.e. original well A1
         - sample 2 would be plate1-A2, i.e. original well A3
         - sample 3 would be plate1-A3, i.e. original well A5
            ...
         - sample 13 would be plate1-B1, i.e. original well C1
         - sample 14 would be plate1-B2, i.e. original well C3
            ...
         - sample 97 would be plate1-H11, i.e. original well O21
         - sample 96 would be plate1-H12, i.e. original well O23
            ...
         - sample 97 would be plate2-A1, i.e. original well A2
         - sample 98 would be plate2-A2, i.e. original well A4
            ...
         - sample 383 would be plate4-H11, i.e. original well P22
         - sample 384 would be plate4-H12, i.e. original well P24
    """
    # first convert everything to 0-based - it makes all the math MUCH simpler
    sample_number -= 1
    # first 96 are plate1, then plate2, plate3, plate4.
    plate_number_96plate = sample_number // 96
    well_number_96plate =  sample_number %  96
    # split well number into row/column number (96-well plates are 8x12)
    row_number_96plate = well_number_96plate // 12
    col_number_96plate = well_number_96plate %  12
    # now get the original 384-well rows/columns, based on the staggering
    if_second_row = plate_number_96plate > 1
    if_second_col = plate_number_96plate % 2
    row_number_384plate = row_number_96plate*2 + if_second_row
    col_number_384plate = col_number_96plate*2 + if_second_col
    assert 0 <= row_number_384plate < 16
    assert 0 <= col_number_384plate < 24
    return '%s%02d'%(ascii_uppercase[row_number_384plate], col_number_384plate+1)


############################################ Conversion from readcounts to codewords ############################################

# MAYBE-TODO should this be in mutant_utilities or something? Not really deconvolution-specific...
def datasets_to_readcount_table(insertion_pool_joint_dataset, dataset_names_in_order, use_perfect_reads=False):
    """ Given an insertional pool multi-dataset, return a insertion_position:readcount_list dict and a list of datasets in order.

    Insertion_pool_joint_dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, multi-dataset, 
     with each member dataset corresponding to a pool.
    dataset_names_in_order should be a list of the dataset names in the order in which they should be listed in the outputs.

    The two outputs are:
     1) a insertion_position:readcount_list, with keys being the positions of all the insertions in the original multi-dataset, 
        and the values being lists of readcounts for that insertion position in each single dataset, 
        in the order given by dataset_names_in_order.

    If use_perfect_reads, the values in the final table will be perfect read counts; otherwise total read counts.
    """
    # MAYBE-TODO implement option for taking a list of datasets rather than a single joint dataset, too?  (would then also have to take a list of names, which I guess is the same as dataset_names_in_order)
    insertion_readcount_table = {}
    for insertion in insertion_pool_joint_dataset:
        if use_perfect_reads:   readcounts = [insertion.by_dataset[dataset].perfect_read_count for dataset in dataset_names_in_order]
        else:                   readcounts = [insertion.by_dataset[dataset].total_read_count for dataset in dataset_names_in_order]
        insertion_readcount_table[insertion.position] = readcounts
    return insertion_readcount_table


# MAYBE-TODO should this be in mutant_utilities or something? Not really deconvolution-specific...
def normalize_readcount_table(insertion_readcount_table, multiplier=10**6):
    """ Given an insertion:dataset_readcount_list dict, normalize the readcounts to the totals for each dataset, IN PLACE.

    Input should be same as the output from datasets_to_readcount_table.  Input will be modified in-place.
    """
    dataset_readcount_totals = [sum(all_readcounts) for all_readcounts in zip(*insertion_readcount_table.values())]
    for insertion_pos, dataset_readcounts in list(insertion_readcount_table.items()):
        insertion_readcount_table[insertion_pos] = [readcount/total * multiplier
                                                   for readcount,total in zip(dataset_readcounts,dataset_readcount_totals)]
    # doesn't return anything - modifies the input in place.
    # TODO unit-test!


def readcounts_to_presence__cutoffs(readcount_list, cutoffs):
    """ Given a readcount list and a cutoff list (or one cutoff), return True/False list for whether readcount>=cutoff.
    """
    # if cutoffs is a single value, make it a list repeating that single value
    try:                len(cutoffs)
    except TypeError:   cutoffs = [cutoffs for _ in range(len(readcount_list))]
    if not len(cutoffs)==len(readcount_list):
        raise DeconvolutionError("readcount_list and cutoffs must be the same length!")
    return [readcount>=cutoff for readcount,cutoff in zip(readcount_list, cutoffs)]


def readcounts_to_presence__mutant_1(readcount_list, N_always_present, N_always_absent, overall_min=2, 
                                     present_level_function=numpy.median, absent_level_function=numpy.median, 
                                     cutoff_position=0.5, min_cutoff=1):
    """ Given a readcount list for a mutant, use a complex method to decide which readcounts indicate presence/absence.

    ___
    """
    # TODO give more details in docstring!
    # TODO should cutoff_positon be a function instead?
    if N_always_present+N_always_absent > len(readcount_list):
        raise DeconvolutionError("N_always_present+N_always_absent should never be more than len(readcount_list)!")
    if N_always_present<1 or N_always_absent<0:
        raise DeconvolutionError("N_always_present must be >0, and N_always_absent >=0!")
    sorted_readcounts = sorted(readcount_list)
    if N_always_absent == 0:
        absent_level = 0
    else:
        absent_readcounts = sorted_readcounts[:N_always_absent]
        absent_level = absent_level_function(absent_readcounts)
    present_readcounts = sorted_readcounts[-N_always_present:]
    present_level = present_level_function(present_readcounts)
    # TODO how exactly should the overall_min be applied?  
    #  to throw away entire codewords as all-0, should we be comparing it to present_level, or to max?
    if present_level < overall_min:
        return [False for readcount in readcount_list]
    if absent_level >= present_level:
        raise DeconvolutionError("Absent level >= present level - should never happen! (%s and %s, based on %s and %s)"%(
            absent_level, present_level, ([] if N_always_absent==0 else absent_readcounts), present_readcounts))
    cutoff = absent_level + cutoff_position*(present_level-absent_level)
    # adjust cutoff
    cutoff = max(cutoff, min_cutoff)
    # MAYBE-TODO what if cutoff ends up below max(absent_readcounts), or above min(present_readcounts)?  Fix?
    # MAYBE-TODO should there be some minimum difference between absent_level and present_level?
    return [readcount>=cutoff for readcount in readcount_list]


def readcounts_to_codewords(insertion_readcount_table, conversion_function):
    """ Given insertion:reads-per-pool data, return insertion:codeword dict showing which insertion was present in which pool. 

    Insertion_readcount_table should be in the same format as the output of datasets_to_readcount_table.
    Conversion_function should be a function that takes a readcount list and returns a True/False list indicating 
     which of those readcounts are considered present/absent (choose one of the readcounts_to_presence__* functions above, 
      and do a lambda to set the other arguments).  

    The output is an insertion_position:presence_codeword dictionary with the same shape as the first input, 
     with presence_codeword a binary_code_utilities.Binary_codeword 
     representing whether a given insertion is considered present in each pool.
    """
    insertion_codeword_dict = {}
    for insertion_pos, readcount_list in insertion_readcount_table.items():
        presence_list = conversion_function(readcount_list)
        insertion_codeword_dict[insertion_pos] = binary_code_utilities.Binary_codeword(presence_list)
    return insertion_codeword_dict


########################################## Matching observed to expected codewords ##############################################

def read_codewords_from_file(infile_name, new_sample_names=None):
    """ Read sample codewords used for combinatorial pooling, from file generated by robotic_plate_transfer.write_data_to_outfile.

    (robotic_plate_transfer is in the ../../combinatorial_pooling/code folder.)
    Only reads the first "sample_number plate_and_well_position codeword transfers volume" table in the infile, ignores the rest.

    Returns sample_number:codeword dictionary, where codewords are binary_code_utilities.Binary_codeword objects.

    If new_sample_names isn't None, it should be an old_number:new_name dictionary - the new names will then be used in the output.
    """
    if new_sample_names is not None:
        if len(set(new_sample_names.values())) != len(new_sample_names):
            raise DeconvolutionError("All values in the new_sample_names dict must be unique!")
    sample_codewords = {}
    inside_relevant_section = False
    for line in open(infile_name):
        # skip lines until inside the section we want to parse
        if not inside_relevant_section:
            if line.startswith('sample_number\tplate_and_well_position\tcodeword'):
                inside_relevant_section=True
                continue
        # inside the section we want to parse:
        if inside_relevant_section:
            # an empty line or another header means the end of the section - just finish the loop.
            if (not line.strip()) or line.startswith('sample_number'):
                break
            # for each line in the section, just grab sample name and codeword from each line and put in the data dict
            fields = line.strip('\n').split('\t')
            try:                sample_name, codeword_string = fields[0], fields[2]
            except IndexError:  raise DeconvolutionError("Cannot parse line! \"%s\""%line)
            # optionally convert original to new sample names
            if new_sample_names is not None:
                sample_name = new_sample_names[sample_name]
            sample_codewords[sample_name] = binary_code_utilities.Binary_codeword(codeword_string)
    return sample_codewords


def find_closest_sample_codeword(ref_codeword, sample_codewords, min_distance_difference=1):
    """ Returns the sample with the closest codeword match to the ref_codeword, and the Hamming distance.

    Inputs: a 0/1 string codeword (a binary_code_utilities.Binary_codeword instance, or just a string), 
    and a sample:codeword dictionary containing more such codewords (values must be unique).

    Finds the sample(s) with codeword(s) with the lowest and second-lowest Hamming distance to the ref_codeword:
        - if there is a single sample with the lowest distance, and the difference between the lowest and second-lowest distances 
            is <=min_distance_difference, return the (closest_sample_name, min_Hamming_distance) tuple for the lowest.
        - otherwise, return a (None, min_Hamming_distance) tuple.
    """
    # useful stuff in binary_code_utilities: Binary_codeword object, Hamming_distance, bit_change_count (separate 0->1 and 1->0)
    # make sure codewords are unique
    if sample_codewords is not None:
        if len(set(sample_codewords.values())) != len(sample_codewords):
            raise DeconvolutionError("All values in the sample_codewords dict must be unique!")
    if min_distance_difference<1:
        raise DeconvolutionError("min_distance_difference must be an integer 1 or higher!")
    # LATER-TODO should probably rename Hamming_distance to be lowercase, since it's a function...
    if not isinstance(ref_codeword, binary_code_utilities.Binary_codeword):
        ref_codeword = binary_code_utilities.Binary_codeword(ref_codeword)
    # Just calculate the Hamming distance to all expected codewords
    #  MAYBE-TODO could probably optimize this a lot!  If only with caching...
    sample_to_distance = {sample:binary_code_utilities.Hamming_distance(ref_codeword, codeword) 
                                   for sample,codeword in sample_codewords.items()}
    min_distance = min(sample_to_distance.values())
    low_dist_samples = [sample for sample,distance in sample_to_distance.items() 
                        if distance < min_distance+min_distance_difference]
    if len(low_dist_samples)==1:    return (low_dist_samples[0], min_distance) 
    else:                           return (None, min_distance) 
    # TODO is this what I actually want to return, or something else?...  Could optionally return the full sorted distance_to_N_samples, or the top 2 of that, or something...


def match_insertions_to_samples(insertion_codewords, sample_codewords, min_distance_difference=1):
    """ Given observed insertion codewords and expected sample codewords, find best sample codeword match for each insertion.

    First two inputs should be insertion_pos:codeword and sample_name:codeword dictionaries, 
     with codewords beint binary_code_utilities.Binary_codeword instances.
    The find_closest_sample_codeword function is used to match insertion to sample codewords; 
     the min_distance_difference arg should be a number - see find_closest_sample_codeword for how it works. 

    The outputs are an insertion_pos:sample_name and an insertion_pos:min_distance dictionary.
    """
    # TODO add option for how many errors should be allowed at all, rather than just always returning the closest unique match, even if it has 4 errors!  Maybe allow splitting that into 0->1 and 1->0 errors separately.
    insertion_samples = {}
    insertion_codeword_distances = {}
    for insertion_pos, insertion_codeword in insertion_codewords.items():
        best_match_sample, min_distance = find_closest_sample_codeword(insertion_codeword, sample_codewords, min_distance_difference)
        insertion_samples[insertion_pos] = best_match_sample
        insertion_codeword_distances[insertion_pos] = min_distance
    return insertion_samples, insertion_codeword_distances
    # TODO unit-test!


############################################### General deconvolution functions ##################################################

def combinatorial_deconvolution(insertion_pool_joint_dataset, sample_codeword_filename, pool_names_in_order, conversion_function, 
                                min_distance_difference=1, new_sample_names=None, 
                                normalize_pool_readcounts=True, normalization_multiplier=10**6, use_perfect_reads=False, 
                                return_readcounts=False):
    """ ___
    """
    # TODO write docstring! Avoid repetition though.
    # TODO wait, what's the relationship between pool_names_in_order and new_sample_names??
    sample_codewords = read_codewords_from_file(sample_codeword_filename, new_sample_names)
    insertion_readcount_table = datasets_to_readcount_table(insertion_pool_joint_dataset, pool_names_in_order, use_perfect_reads)
    if normalize_pool_readcounts: normalize_readcount_table(insertion_readcount_table, normalization_multiplier)
    insertion_codewords = readcounts_to_codewords(insertion_readcount_table, conversion_function)
    insertion_samples, insertion_codeword_distances = match_insertions_to_samples(insertion_codewords, sample_codewords, 
                                                                                  min_distance_difference)
    # MAYBE-TODO add more options for what is and isn't returned?  Or are those a bad idea?
    return_data = (insertion_codewords, insertion_samples, insertion_codeword_distances)
    if return_readcounts:   return_data = (insertion_readcount_table,) + return_data
    return return_data
    # TODO unit-test!


def get_deconvolution_summary(insertion_samples, insertion_codeword_distances):
    """ Given the outputs of match_insertions_to_samples or combinatorial_deconvolution, generate some summary numbers.

    Inputs are insertion_pos:sample_name and insertion_pos:min_codeword_distance dictionaries; 
        sample_name should be None if multiple samples matched.

    Outputs are the total numbers of matched and unmatched insertions, 
     and min_distance:N_insertions dicts for matched and unmatched insertions.
    """
    matched_insertion_distances = {i:d for i,d in insertion_codeword_distances.items() if insertion_samples[i] is not None}
    unmatched_insertion_distances = {i:d for i,d in insertion_codeword_distances.items() if insertion_samples[i] is None}
    N_total_insertions = len(insertion_samples)
    N_matched_insertions = len(matched_insertion_distances)
    N_unmatched_insertions = len(unmatched_insertion_distances)
    matched_insertion_counts_by_distance = collections.Counter(matched_insertion_distances.values())
    unmatched_insertion_counts_by_distance = collections.Counter(unmatched_insertion_distances.values())
    return N_matched_insertions, N_unmatched_insertions, matched_insertion_counts_by_distance, unmatched_insertion_counts_by_distance
    # TODO unit-test


def print_deconvolution_summary(description, N_matched_insertions, N_unmatched_insertions, matched_insertion_counts_by_distance, 
                                unmatched_insertion_counts_by_distance):
    """ Nicely print the data generated by get_deconvolution_summary.
    """
    # MAYBE-TODO add a collapse_higher_than arg, so that "All insertion counts with min distances >collapse_higher_than will be counted together instead of separately."? (put that bit in docstring)
    total_insertions = N_matched_insertions + N_unmatched_insertions
    print "%s: %s total insertions, %s uniquely matched to a sample, %s matched to 2+ samples (unmatched)."%(description, 
                total_insertions, general_utilities.value_and_percentages(N_matched_insertions, [total_insertions]), 
                general_utilities.value_and_percentages(N_unmatched_insertions, [total_insertions]))
    print " * matched insertions by codeword distance (% of matched, % of all): "
    print "    %s"%(', '.join(['%s: %s'%(d,general_utilities.value_and_percentages(n, [N_matched_insertions, total_insertions]))
                                         for d,n in sorted(matched_insertion_counts_by_distance.items())]))
    print " * unmatched insertions by codeword distance (% of unmatched, % of all): "
    print "    %s"%(', '.join(['%s: %s'%(d,general_utilities.value_and_percentages(n, [N_unmatched_insertions, total_insertions]))
                                         for d,n in sorted(unmatched_insertion_counts_by_distance.items())]))


def pick_best_parameters(deconv_data_by_parameters, N_errors_allowed, N_good_factor, percent_good_factor, max_N_high=None):
    """ Pick the best deconvolution parameters, depending on what aspect of the result we want to optimize.

    The deconv_data_by_parameters arg should be a dict with (N_high, N_low, cutoff_pos, overall_min) tuples as keys and 
     (insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances, summary) tuples as values.
        (The first four should be outputs from combinatorial_deconvolution with return_readcounts=True, or the equivalent;
         summary should be the output from get_deconvolution_summary.)

    Look at the parameters/results in deconv_data_by_parameters.  
    Pick the best set of parameters, optimizing a combination of two factors:
        1) the total number of insertion positions mapped to expected codewords with N_errors_allowed or fewer errors
        2) the percentage of those out of all uniquely mapped insertion positions.
    These two numbers will be multiplied by N_good_factor and percent_good_factor respectively, 
     and the resulting combination will be  maximized.

    The output will be the deconv_data_by_parameters key that gives the optimal result.
    """
    parameters_to_result = {}
    for parameters, deconv_data in deconv_data_by_parameters.items():
        N_matched, N_unmatched, N_matched_by_dist, _ = summary = deconv_data[-1]
        N_good_matched = sum(N_matched_by_dist[x] for x in range(N_errors_allowed+1))
        percent_good_matched = N_good_matched/N_matched*100
        final_result = N_good_matched*N_good_factor + percent_good_matched*percent_good_factor
        if max_N_high is None or parameters[0]<=max_N_high:
            parameters_to_result[parameters] = final_result
    best_parameters = max(parameters_to_result.items(), key=lambda (p,r): r)[0]
    return best_parameters
        

# TODO write function to print deconvolution data!  Fields: plate, well, best/good/decent category for plate and well (if there's more than one), N_errors for plate and well, maybe average_readcount for plate and well, insertion position, probably gene/feature/annotation etc.  Also need to output insertions that were mapped to a plate but not a well, and vice versa - separate files?  Or just give them "-" in place of plate/well in the same file?



################################################## Plotting the data #######################################################

def absent_present_readcounts(readcounts, codeword):
    """ Given a list of readcounts and a binary codeword, return separate lists of readcounts with 0s and with 1s in the codeword.
    """
    readcounts_and_presence = zip(readcounts, codeword.list())
    absent_readcounts = [readcount for (readcount,presence) in readcounts_and_presence if not presence]
    present_readcounts = [readcount for (readcount,presence) in readcounts_and_presence if presence]
    return absent_readcounts, present_readcounts

def plot_mutants_and_cutoffs(insertion_readcount_table, insertion_codewords, order='by_average_readcount', min_readcount=20, 
                             N_rows=60, N_columns=13, N_insertions=None, x_scale='symlog', markersize=4, marker='d', 
                             filename=None, title='insertion readcounts (white - absent, black - present)'):
    """ ____

    If filename is not given, plots will be interactive; otherwise they won't, and will be saved as multiple files.
    """
    # TODO write docstring!
    # MAYBE-TODO add option to only plot insertions present in a specific dataset/list?
    if filename is not None:
        mplt.ioff()
    # TODO add give option for insertion_codewords to be None - all readcounts treated the same, no present/absent distinction
    # make a joint [(readcounts, codeword)] list, discarding insertion position data
    # originally the insertion order is by position - reorder it if desired
    # MAYBE-TODO add different order options?
    readcounts_and_codewords = [(readcounts, insertion_codewords[insertion_pos])
                                 for (insertion_pos,readcounts) in sorted(insertion_readcount_table.items())]
    if order=='by_position':            pass
    elif order=='random':               random.shuffle(readcounts_and_codewords)
    elif order=='by_average_readcount':  readcounts_and_codewords.sort(key=lambda (r,c): numpy.mean(r))
    elif order=='by_median_readcount':  readcounts_and_codewords.sort(key=lambda (r,c): numpy.median(r))
    else:                               raise DeconvolutionError("Unknown order value %s!"%order)
    # if desired, remove any data points that have no readcounts above min_readcount AND have empty codewords
    if min_readcount is not None:       
        readcounts_and_codewords = filter(lambda (r,c): max(r)>=min_readcount or c.weight(), 
                                          readcounts_and_codewords)
    # if desired, take some number of insertions, from the start
    if N_insertions is not None:        readcounts_and_codewords = readcounts_and_codewords[:N_insertions]
    N_per_page = N_rows*N_columns
    N_pages = int(math.ceil(len(readcounts_and_codewords)/N_per_page))
    max_readcount = max(max(readcounts) for (readcounts,_) in readcounts_and_codewords)
    dot_plot_kwargs = dict(linestyle='None', markeredgecolor='k', markersize=markersize)
    for page_N in range(N_pages):
        mplt.figure(figsize=(16,12))
        mplt.suptitle(title + '; page %s/%s'%(page_N+1, N_pages), y=0.92)      # y position to avoid wasting space
        page_first_ins_N = page_N*N_per_page
        for column_N in range(N_columns):
            col_first_ins_N = page_first_ins_N + column_N*N_rows
            if col_first_ins_N >= len(readcounts_and_codewords):
                break
            mplt.subplot(1, N_columns, column_N+1)
            # plot a grey line for each insertion, with dots showing readcounts (filled if "present", unfilled if "absent")
            # TODO currently we don't explicitly know the cutoffs, just the codewords!  Should I calculate (or save) the real cutoffs used, instead?  Or infer them based on the highest absent and lowest present readcount and plot those?
            mplt.hlines(range(min(N_rows, len(readcounts_and_codewords)-col_first_ins_N)), 0, max_readcount, colors='0.7')
            # accumulate all dot positions before plotting them all together for the whole column
            absent_dot_positions_x, absent_dot_positions_y = [], []
            present_dot_positions_x, present_dot_positions_y = [], []
            for row_N in range(N_rows):
                curr_ins_N = col_first_ins_N+row_N
                if curr_ins_N >= len(readcounts_and_codewords):
                    break
                absent_readcounts, present_readcounts = absent_present_readcounts(*readcounts_and_codewords[curr_ins_N])
                absent_dot_positions_x.extend(absent_readcounts)
                absent_dot_positions_y.extend([row_N for _ in absent_readcounts])
                present_dot_positions_x.extend(present_readcounts)
                present_dot_positions_y.extend([row_N for _ in present_readcounts])
            mplt.plot(absent_dot_positions_x, absent_dot_positions_y, marker, markerfacecolor='w', **dot_plot_kwargs)
            mplt.plot(present_dot_positions_x, present_dot_positions_y, marker, markerfacecolor='k', **dot_plot_kwargs)
            mplt.xscale(x_scale)
            # remove the axes/etc we don't want, make things pretty
            mplt.xlim(-1, max_readcount*2)
            mplt.ylim(-1, N_rows)
            mplt.yticks([], [])
            mplt.xticks(mplt.xticks()[0], [])
            ax = mplt.gca()
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.yaxis.set_ticks_position('none')
            mplt.draw()
            # TODO tighten layout
            mplt.subplots_adjust(wspace=0.1)
            #mplt.tight_layout() - this isn't great, messes up suptitle
        if filename is not None:
            basename, ext = os.path.splitext(filename)
            pagenum_len = len(str(N_pages+1))
            format_str = '_%%0%sd'%pagenum_len
            plotting_utilities.savefig(basename + format_str%(page_N+1) + ext)
            mplt.close()
    if filename is not None:
        mplt.ion()


def plot_parameter_space_1(deconv_data_by_parameters, N_errors_allowed, chosen_parameters, info='', 
                           marker='x', alpha=1, colors_by_N_high='orangered black blue'.split()):
    """ Given deconvolution data for various parameter sets, plot # and % of matched with <= N errors, and highlight chosen set.

    The deconv_data_by_parameters arg should be a dict with (N_high, N_low, cutoff_pos, overall_min) tuples as keys and 
     (insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances, summary) tuples as values.
        (The first four should be outputs from combinatorial_deconvolution with return_readcounts=True, or the equivalent;
         summary should be the output from get_deconvolution_summary.)
        
    Each item of deconv_data_by_parameters will be plotted as a dot. 
    The chosen_parameters should be one of the deconv_data_by_parameters keys - that dot will be highlighted.

    N_errors_allowed should be an integer 0 or higher, and will be used to define what the axes are:
        - x axis will the the number of insertion positions that were mapped to expected codewords with N or fewer errors
        - y axis will be that number divided by the total number of uniquely mapped insertion positions (with any #errors)

    Info will be used in the plot title.

    The marker, alpha, and colors_by_N_high args set the physical appearance of the markers.

    """
    if not chosen_parameters in deconv_data_by_parameters.keys():
        raise DeconvolutionError("Chosen parameters not present in result set! %s, %s"%(chosen_parameters, 
                                                                                        deconv_data_by_parameters.keys()))
    N_good_matched_val_dict, percent_good_matched_val_dict = defaultdict(list), defaultdict(list)
    for parameters, deconv_data in deconv_data_by_parameters.items():
        N_matched, N_unmatched, N_matched_by_dist, _ = summary = deconv_data[-1]
        N_good_matched = sum(N_matched_by_dist[x] for x in range(N_errors_allowed+1))
        percent_good_matched = N_good_matched/N_matched*100
        N_good_matched_val_dict[parameters[0]].append(N_good_matched)
        percent_good_matched_val_dict[parameters[0]].append(percent_good_matched)
        if parameters==chosen_parameters:
            chosen_N_good_matched, chosen_percent_good_matched = N_good_matched, percent_good_matched
    min_N_high = min(N_good_matched_val_dict.keys())
    for N_high,N_good_matched_vals in N_good_matched_val_dict.items():
        mplt.plot(N_good_matched_vals, percent_good_matched_val_dict[N_high], 
                  marker=marker, color=colors_by_N_high[N_high-min_N_high], markersize=5, alpha=alpha, linestyle='None')
    mplt.plot(chosen_N_good_matched, chosen_percent_good_matched,  
              marker='o', markerfacecolor='None', markeredgecolor='r', markersize=10, markeredgewidth=2, linestyle='None')
    mplt.title("Deconvolution results for different parameters%s.\n Allowing %s errors; red circle marks chosen parameters."%(
        ': %s'%info if info else '', N_errors_allowed))
    mplt.xlabel("number of insertions mapped with at most %s errors"%N_errors_allowed)
    mplt.ylabel("% of those out of all uniquely mapped insertions")
    # MAYBE-TODO could try changing the marker shape/size/color to reflect the other parameters?


###################################################### Testing ###########################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__get_original_384_well_numbers(self):
        # assuming the 384-well plate was split into four 96-well plates, every other row/column - SEE DOCSTRING FOR INFO
        #   all the examples here are straight from the get_original_384_well_numbers docstring.
        assert get_original_384_well_numbers(1)   == 'A01'
        assert get_original_384_well_numbers(2)   == 'A03'
        assert get_original_384_well_numbers(3)   == 'A05'
        assert get_original_384_well_numbers(13)  == 'C01'
        assert get_original_384_well_numbers(14)  == 'C03'
        assert get_original_384_well_numbers(95)  == 'O21'
        assert get_original_384_well_numbers(96)  == 'O23'
        assert get_original_384_well_numbers(97)  == 'A02'
        assert get_original_384_well_numbers(98)  == 'A04'
        assert get_original_384_well_numbers(383) == 'P22'
        assert get_original_384_well_numbers(384) == 'P24'

    def _make_positions(self, N=3):
        """ Help function - make three insertion positions, on chromosomes 1, 2, 3, varying strands/etc, immutable. """
        positions = []
        _I = mutant_analysis_classes.Insertion_position
        positions.append(_I('chr1', '+', position_before=100, immutable=True))
        positions.append(_I('chr2', '-', position_after=501, immutable=True))
        positions.append(_I('chr3', '+', position_before=200, position_after=201, immutable=True))
        positions.append(_I('chr4', '-', position_after=20, immutable=True))
        positions.append(_I('chr5', '_', position_before=444, immutable=True))
        return positions[:N]

    def test__datasets_to_readcount_table(self):
        # make a test case with 3 pools and 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        pools = ['A', 'B', 'C']
        pos1, pos2, pos3 = self._make_positions()
        insertion1 = mutant_analysis_classes.Insertional_mutant(pos1, multi_dataset=True)
        insertion2 = mutant_analysis_classes.Insertional_mutant(pos2, multi_dataset=True)
        insertion3 = mutant_analysis_classes.Insertional_mutant(pos3, multi_dataset=True)
        insertions = [insertion1, insertion2, insertion3]
        # the three numerical arguments to add_counts are total_reads, perfect_reads, sequence_variants (not used).
        insertion1.add_counts(0,  0,  0, dataset_name=pools[0])
        insertion1.add_counts(1,  0,  1, dataset_name=pools[1])
        insertion1.add_counts(10, 9,  1, dataset_name=pools[2])
        insertion2.add_counts(9,  7,  1, dataset_name=pools[0])
        insertion2.add_counts(20, 20, 1, dataset_name=pools[1])
        insertion2.add_counts(0,  0,  0, dataset_name=pools[2])
        insertion3.add_counts(1,  0,  1, dataset_name=pools[0])
        insertion3.add_counts(0,  0,  0, dataset_name=pools[1])
        insertion3.add_counts(3,  3,  1, dataset_name=pools[2])
        dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset=True)
        for insertion in insertions:  dataset.add_mutant(insertion)
        assert datasets_to_readcount_table(dataset, pools, use_perfect_reads=False) == {pos1:[0,1,10], pos2:[9,20,0], pos3:[1,0,3]}
        assert datasets_to_readcount_table(dataset, pools, use_perfect_reads=True) ==  {pos1:[0,0,9],  pos2:[7,20,0], pos3:[0,0,3]}

    def _make_basic_readcount_table(self):
        # make a test case with 3 pools and 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        pos1, pos2, pos3 = self._make_positions()
        return {pos1:[0,1,10], pos2:[9,20,0], pos3:[1,0,3]}

    def _simplify_codeword_output(self, insertion_codeword_dict, ins_positions):
        """ Convenience function for easier processing of codeword tables - return them as a string in given order. """
        return ' '.join([str(insertion_codeword_dict[ins_position]) for ins_position in ins_positions])

    # Note: all the readcounts_to_presence__* tests below ALSO test readcounts_to_codewords, 
    #  since that's just a convenience function to do readcounts_to_presence on a larger scale and change outputs to codewords.

    def test__readcounts_to_presence__cutoffs(self):
        # 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        insertions = self._make_positions()
        readcounts = self._make_basic_readcount_table()
        _S = self._simplify_codeword_output
        function_with_cutoff = lambda cutoff: (lambda x: readcounts_to_presence__cutoffs(x, cutoff))
        # single cutoff (checking all relevant values)
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(0)), insertions)  == '111 111 111'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(1)), insertions)  == '011 110 101'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(2)), insertions)  == '001 110 001'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(3)), insertions)  == '001 110 001'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(4)), insertions)  == '001 110 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(9)), insertions)  == '001 110 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(10)), insertions) == '001 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(11)), insertions) == '000 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(20)), insertions) == '000 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(21)), insertions) == '000 000 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(99)), insertions) == '000 000 000'
        # per-pool cutoffs
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([1,1, 1])), insertions) == '011 110 101'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([3,10,1])), insertions) == '001 110 001' 
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([2,10,4])), insertions) == '001 110 000' 

    def test__readcounts_to_presence__mutant_1(self):
        # 3 insertions with various readcount numbers: perfect, noisy, absent, weird
        # MAYBE-TODO try changing the order of the readcounts and making sure the results are the same (but in different order)?
        insertions = self._make_positions(4)
        readcounts = dict(zip(insertions, ([0,0,0,0,10,10,10], [0,1,1,3,8,10,15], [0,0,0,0,1,1,1], [0,1,2,3,4,5,6])))
        _S = self._simplify_codeword_output
        # arguments to readcounts_to_presence__mutant_1 are: readcount_list, N_always_present, N_always_absent, overall_min, 
        #                                               present_level_function, absent_level_function, cutoff_position, min_cutoff
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 1, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000111 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 2, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # upping the overall_min makes more all-0 results
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 6, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000000'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 10, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000000'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 11, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000000 0000000 0000000 0000000'
        # upping min_cutoff makes more single 0s
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 5)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000011'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 10)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000011 0000000 0000000'
        # changing cutoff_position up/down from 0.5 - more/fewer 0s
        #  (readcounts from above: [0,0,0,0,10,10,10], [0,1,1,3,8,10,15], [0,0,0,0,1,1,1], [0,1,2,3,4,5,6])
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.8, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000011 0000000 0000011'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.2, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0001111 0000000 0011111'
        ### changing the level functions from median to max/min:
        # mean of 2 numbers is the same as median; not true for 3 numbers, but in this case it ends up the same
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 2, numpy.mean, numpy.mean, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.mean, numpy.mean, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # min/max - actually ends up pretty boring...
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, max, min, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, min, max, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # TODO make more interesting min/max tests!  Or other functions
        # TODO test N_always_absent == 0

    def test__read_codewords_from_file(self):
        # convenience function: compare real output (includes Binary_codeword objects) to simple string representation of dict
        def _compare(output, expected_string):
            assert len(output) == len(expected_string.split(' '))
            for single_string in expected_string.split(' '):
                key,val = single_string.split(':')
                assert str(output[key]) == val
        # basic file
        infile1 = 'test_data/INPUT_codewords_1.txt'
        output1 = read_codewords_from_file(infile1)
        _compare(output1, '0:001 1:010 2:011 3:100 4:101 5:110 6:111')
        # a file with another code, and with a mirror-codeword section too
        infile2 = 'test_data/INPUT_codewords_2.txt'
        output2 = read_codewords_from_file(infile2)
        _compare(output2, '0:0011 1:0101 2:0110 3:1001 4:1010 5:1100 6:1111')
        # trying out the new_sample_names optional dict
        output1b = read_codewords_from_file(infile1, {str(x):str(x+10) for x in range(7)})
        _compare(output1b, '10:001 11:010 12:011 13:100 14:101 15:110 16:111')
        output2b = read_codewords_from_file(infile2, {'0':'dA', '1':'dB', '2':'dC', '3':'dD', '4':'dE', '5':'dF', '6':'dG'})
        _compare(output2b, 'dA:0011 dB:0101 dC:0110 dD:1001 dE:1010 dF:1100 dG:1111')
        # fail if new_sample_names values are non-unique
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):'1' for x in range(7)})
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):min(x,5) for x in range(7)})
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):('A' if x<3 else 'B') for x in range(7)})

    def test__find_closest_sample_codeword(self):
        sample_codewords = {x[0]:binary_code_utilities.Binary_codeword(x[2:]) for x in 'A:1100 B:1010 C:0011 D:1000'.split()}
        # diff 1 - if the top and second Hamming distance differ by 1, take the top one
        inputs_outputs_diff1 = ('1100 A 0, 1010 B 0, 0011 C 0, 1000 D 0, '
                                +'1111 None 2, 0000 D 1, 0110 None 2, 0111 C 1, 1110 None 1, 0001 C 1, 0100 A 1')
        # diff 1 - if the top and second Hamming distance differ by 1, count that as None - they must differ by at least 2
        inputs_outputs_diff2 = ('1100 None 0, 1010 None 0, 0011 C 0, 1000 None 0, '
                                +'1111 None 2, 0000 None 1, 0110 None 2, 0111 C 1, 1110 None 1, 0001 None 1, 0100 None 1')
        for diff, inputs_outputs in [(1, inputs_outputs_diff1), (2, inputs_outputs_diff2)]:
            for input_output_str in inputs_outputs.split(', '):
                input_str, sample, distance = input_output_str.split(' ')
                distance = int(distance)
                for input_val in (input_str, binary_code_utilities.Binary_codeword(input_str)):
                    out_sample, out_dist = find_closest_sample_codeword(input_val, sample_codewords, diff)
                    if sample == 'None':    assert (out_sample is None and out_dist == distance)
                    else:                   assert (out_sample, out_dist) == (sample, distance)

    # LATER-TODO add more unit-tests!



if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
