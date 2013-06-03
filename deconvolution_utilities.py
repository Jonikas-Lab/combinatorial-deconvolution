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
import collections
import unittest
# other packages
# my modules
import general_utilities
import binary_code_utilities
import mutant_analysis_classes

class DeconvolutionError(Exception):
    """ Exceptions in the deconvolution_utilities module."""
    pass


def readcounts_to_presence(insertion_pool_joint_dataset, one_cutoff=None, cutoff_per_pool=None, cutoff_per_insertion=None):
    """ Given a reads-per-insertion dataset and cutoffs, determine which insertion was above/below cutoff in which pool. 

    Insertion_pool_joint_dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, multi-dataset, 
     with each member dataset corresponding to a pool.

    There are three ways of specifying readcount cutoffs, and exactly one should be provided:
        - one_cutoff: a single number that will be used for all pool/insertion combinations
        - cutoff_per_pool: a pool_name:number dict, so the cutoff depends on the pool
        - cutoff_per_insertion: a insertion_position:number dict, so the cutoff depends on the insertion
    
    The output is an insertion_position:pool_name:X nested dictionary, with X 0 if the insertion had <cutoff reads in the pool, 
     and 1 if it had >=cutoff reads.
    """
    if sum(x is not None for x in (one_cutoff, cutoff_per_pool, cutoff_per_insertion)) != 1:
        raise DeconvolutionError("Must provide exactly one of the cutoff arguments!")
    if one_cutoff is not None:           _get_cutoff = lambda pool_name, insertion_position: one_cutoff
    if cutoff_per_pool is not None:      _get_cutoff = lambda pool_name, insertion_position: cutoff_per_pool[pool_name]
    if cutoff_per_insertion is not None: _get_cutoff = lambda pool_name, insertion_position: cutoff_per_insertion[insertion_position]
    # TODO there should also be an option with the cutoff depending on BOTH pool and insertion, somehow... Normalize the per-pool readcounts before applying the per-insertion cutoff, or something?
    insertion_pool_presence_dict = {}
    for insertion in insertion_pool_joint_dataset:
        insertion_pool_presence_dict[insertion.position] = collections.defaultdict(int)
        for pool_name,pool_insertion_data in insertion.by_dataset.items():
            if_present = int(bool(pool_insertion_data.total_read_count >= _get_cutoff(pool_name, insertion.position)))
            insertion_pool_presence_dict[insertion.position][pool_name] = if_present
    return insertion_pool_presence_dict


def presence_to_codewords(insertion_pool_presence_dict, pool_names_in_order=None, quiet=False):
    """ Convert an insertion/pool absence/presence dict to insertion codeword dict.
    
    Given an insertion_pos:pool_name:0/1 dict and a list of pool names in order, 
     return insertion_pos:pool_codeword dict, where pool_codeword is a binary_code_utilities.Binary_codeword object
      based on the string of 0/1 values for the given insertion, for pools as given in pool_names_in_order.

    pool_names_in_order does NOT have to include all the pool names in insertion_pool_presence_dict, 
     but it cannot include any pools that aren't in it.
    If pool_names_in_order is None, make it simply the sorted list of the pool names from insertion_pool_presence_dict.
    """
    # if pool order not given, take the basic sorted order (assuming each insertion_pool_presence_dict value has same keys)
    if pool_names_in_order is None:
        all_pool_names = set.union(*[set(pool_presence_dict.keys()) for pool_presence_dict in insertion_pool_presence_dict.values()])
        pool_names_in_order = sorted(all_pool_names)
        if not quiet:
            print "Inferring pool name order: %s"%(', '.join(pool_names_in_order))
    # convert the 0/1 values to codewords.
    insertion_codewords = {}
    for insertion_pos, pool_presence_dict in insertion_pool_presence_dict.items():
        codeword_string = ''.join(str(pool_presence_dict[pool_name]) for pool_name in pool_names_in_order)
        insertion_codewords[insertion_pos] = binary_code_utilities.Binary_codeword(codeword_string)
    return insertion_codewords
    # TODO unit-test!


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
    insertion_samples = {}
    insertion_codeword_distances = {}
    for insertion_pos, insertion_codeword in insertion_codewords.items():
        best_match_sample, min_distance = find_closest_sample_codeword(insertion_codeword, sample_codewords, min_distance_difference)
        insertion_samples[insertion_pos] = best_match_sample
        insertion_codeword_distances[insertion_pos] = min_distance
    return insertion_samples, insertion_codeword_distances
    # TODO unit-test!


def combinatorial_deconvolution(insertion_pool_joint_dataset, sample_codeword_filename, 
                                one_cutoff=None, cutoff_per_pool=None, cutoff_per_insertion=None, 
                                pool_names_in_order=None, min_distance_difference=1, new_sample_names=None):
    """ ___
    """
    # TODO write docstring! Avoid repetition though.
    sample_codewords = read_codewords_from_file(sample_codeword_filename, new_sample_names)
    insertion_presence = readcounts_to_presence(insertion_pool_joint_dataset, one_cutoff, cutoff_per_pool, cutoff_per_insertion)
    # MAYBE-TODO do we need to give pool_names_in_order to presence_to_codewords, or is the default all right?
    insertion_codewords = presence_to_codewords(insertion_presence, pool_names_in_order)
    insertion_samples, insertion_codeword_distances = match_insertions_to_samples(insertion_codewords, sample_codewords, 
                                                                                  min_distance_difference)
    return insertion_codewords, insertion_samples, insertion_codeword_distances
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
    print " * matched insertions by codeword distance: %s"%(', '.join(['%s: %s'%(d,n) for d,n 
                                                                       in sorted(matched_insertion_counts_by_distance.items())]))
    print " * unmatched insertions by codeword distance: %s"%(', '.join(['%s: %s'%(d,n) for d,n 
                                                                       in sorted(unmatched_insertion_counts_by_distance.items())]))


###################################################### Testing ###########################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__readcounts_to_presence(self):
        # convenience function for easier output checking:
        def _convert_output(output, pools, insertions):
            """ Given insertion:pool:0/1 output dict and lists of pools and insertions in order, return string like '01 11'. """
            insertion_codewords = {insertion.position: ''.join(str(output[insertion.position][pool]) for pool in pools) 
                                for insertion in insertions}
            return ' '.join([insertion_codewords[insertion.position] for insertion in insertions])
        # make a test case with 3 pools and 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        pools = ['A', 'B', 'C']
        pos1 = mutant_analysis_classes.Insertion_position('chr1', '+', position_before=100, immutable=True)
        pos2 = mutant_analysis_classes.Insertion_position('chr2', '-', position_after=501, immutable=True)
        pos3 = mutant_analysis_classes.Insertion_position('chr3', '+', position_before=200, position_after=201, immutable=True)
        insertion1 = mutant_analysis_classes.Insertional_mutant(pos1, multi_dataset=True)
        insertion2 = mutant_analysis_classes.Insertional_mutant(pos2, multi_dataset=True)
        insertion3 = mutant_analysis_classes.Insertional_mutant(pos3, multi_dataset=True)
        insertions = [insertion1, insertion2, insertion3]
        # the three numerical arguments to add_counts are total_reads,perfect_reads,sequence_variants - only the first matters.
        insertion1.add_counts(0, 0, 0, dataset_name=pools[0])
        insertion1.add_counts(1, 1, 1, dataset_name=pools[1])
        insertion1.add_counts(10, 10, 1, dataset_name=pools[2])
        insertion2.add_counts(9, 9, 1, dataset_name=pools[0])
        insertion2.add_counts(20, 20, 1, dataset_name=pools[1])
        insertion2.add_counts(0, 0, 0, dataset_name=pools[2])
        insertion3.add_counts(1, 1, 1, dataset_name=pools[0])
        insertion3.add_counts(0, 0, 0, dataset_name=pools[1])
        insertion3.add_counts(3, 3, 1, dataset_name=pools[2])
        dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset=True)
        for insertion in insertions:  dataset.add_mutant(insertion)
        # single cutoff (checking all relevant values)
        #   3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=0), pools, insertions)  == '111 111 111'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=1), pools, insertions)  == '011 110 101'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=2), pools, insertions)  == '001 110 001'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=3), pools, insertions)  == '001 110 001'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=4), pools, insertions)  == '001 110 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=9), pools, insertions)  == '001 110 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=10), pools, insertions) == '001 010 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=11), pools, insertions) == '000 010 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=20), pools, insertions) == '000 010 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=21), pools, insertions) == '000 000 000'
        assert _convert_output(readcounts_to_presence(dataset, one_cutoff=99), pools, insertions) == '000 000 000'
        # per-pool cutoffs
        #   3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        _output = lambda CpD: _convert_output(readcounts_to_presence(dataset, cutoff_per_pool=CpD), pools, insertions)
        assert _output({'A':1, 'B':1, 'C':1}) == '011 110 101'
        assert _output({'A':3, 'B':10, 'C':1}) == '001 110 001'
        assert _output({'A':2, 'B':10, 'C':4}) == '001 110 000'
        # per-insertion cutoffs
        #   3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        _output = lambda CpM: _convert_output(readcounts_to_presence(dataset, cutoff_per_insertion=CpM), pools, insertions)
        assert _output({pos1:1, pos2:1, pos3:1}) ==  '011 110 101'
        assert _output({pos1:2, pos2:2, pos3:1}) ==  '001 110 101'
        assert _output({pos1:2, pos2:2, pos3:4}) ==  '001 110 000'
        assert _output({pos1:2, pos2:10, pos3:4}) == '001 010 000'
        # make sure that it only works with exactly one codeword argument - won't work with 0, any combination of 2, or all 3.
        O, D, M = 1, {'A':1, 'B':1, 'C':1}, {pos1:1, pos2:1, pos3:1}
        self.assertRaises(DeconvolutionError, readcounts_to_presence, dataset)
        self.assertRaises(DeconvolutionError, readcounts_to_presence, dataset, one_cutoff=O, cutoff_per_pool=D)
        self.assertRaises(DeconvolutionError, readcounts_to_presence, dataset, one_cutoff=O, cutoff_per_insertion=M)
        self.assertRaises(DeconvolutionError, readcounts_to_presence, dataset, cutoff_per_pool=D, cutoff_per_insertion=M)
        self.assertRaises(DeconvolutionError, readcounts_to_presence,dataset, one_cutoff=O,cutoff_per_pool=D,cutoff_per_insertion=M)

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
