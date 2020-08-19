from simlib import fast_similarity
from Levenshtein import distance

"""Note that the default cost in simlib is 0.55 (as of 08/06/2020); was 0.51"""

max_possible_distance = 100 # this should be the max CDR3 length (for safety); anything greater is overkill
distance_to_similarity_list = [fast_similarity(x) for x in range(max_possible_distance)]

def get_similarity_matrix(unique_seqs_from_all_repertoires):
	
	"""Argument: an iterable (list/tuple) of unique species
	Return: all-against-all similarity matrix based on the fast_similarity method as used in Arora et al. 2018 BiorXiv"""
	similarity_list = ([ distance_to_similarity_list[distance(seq_1, seq_2)] for seq_2 in unique_seqs_from_all_repertoires ] for seq_1 in unique_seqs_from_all_repertoires)
	return similarity_list
