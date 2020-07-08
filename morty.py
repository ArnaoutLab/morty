
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from ast import literal_eval
from collections import Counter, defaultdict
from functools import reduce
from itertools import combinations, product, chain
from Levenshtein import distance
from math import log10, log, ceil, exp
from numpy import dot, array, sum as np_sum, concatenate, identity, savetxt
from numpy.random import choice as np_choice
from operator import itemgetter, mul
from os.path import expanduser, isfile, join
from random import random, choice, shuffle, randint
from simlib import stochastic_similarity, set_cython_seed, fast_similarity
from string import *
from sys import argv, exit
from time import time, sleep, strftime

import numpy as np
import os
import subprocess

""" FUNCTIONS """

def eval_right_side(txt):
	return literal_eval(txt.split("=")[1].str.strip())


def prod(iterable):
	"""The equivalent of sum() for computing the product of the terms
	Arguments:
		iterable: takes an iterable [such as a list] and returns the product of the values"""
	return reduce(mul, iterable, 1)


def old_similarity(average_length, x):
	# average_length = (len(seq_1)+len(seq_2))/2.
	raised_expression = (4./average_length) * x
	old_similarity_measure = 0.1**raised_expression
	return old_similarity_measure


def get_distance_to_similarity_list(cost, method):
	max_possible_distance = 100 # this should be the max CDR3 length (for safety); anything greater is overkill
	if method == "fast_similarity":
		distance_to_similarity_list = [fast_similarity(x, cost) for x in range(max_possible_distance)]

	elif method == "old_similarity":
		# Note below distance_to_similarity_list is a 2D matrix: distance_to_similarity_list[distance][average_length]
		# could also be made a dictionary; an array is probably fastest but this is a list of lists
		# we don't think we'll be using this often
		average_lengths = arange(0, max_possible_distance, 0.5)
		distance_to_similarity_list = []
		for x in range(max_possible_distance):
			dummy_ = []
			for average_length in average_lengths:
				dummy.append(old_similarity(x, average_length))
				distance_to_similarity_list.append(dummy)
	return distance_to_similarity_list
	

# def calculate_similarity(seq_1, seq_2, distance_to_similarity_list=distance_to_similarity_list, method="fast_similarity"):
def calculate_similarity(seq_1, seq_2, cost, method, distance_to_similarity_list):
	
	"""type : fast_similarity vs. old_similarity"""

	# distance_to_similarity_list = get_distance_to_similarity_list(cost, method)
	
	if method=="fast_similarity":
		return distance_to_similarity_list[distance(seq_1, seq_2)]
	
	elif method=="old_similarity":
		x = distance(seq_1, seq_2) # edit distance
		average_length = (len(seq_1)+len(seq_2))/2.
		return distance_to_similarity_list[x][average_length]



def calculate_similarity_BAK(seq_1, seq_2, method="fast_similarity"):

	"""type : fast_similarity vs. old_similarity"""
	
	x = distance(seq_1, seq_2) # edit distance
	
	if method == "fast_similarity": # see log Harry Burke log.rtfd entry 080518RO
		return fast_similarity(x, cost) # Normally this function would accept arguments x and kappa but in simlib, kappa=0.51 by default.
	
	elif method == "old_similarity":
		average_length = (len(seq_1)+len(seq_2))/2.
		return old_similarity(x)


def calculate_alpha_diversity(unique_seqs_from_this_repertoires, p, repertoire_to_seq_prob_hash, alpha_diversity_repertoire, recon_file_for_this_repertoire, list_of_qs, cost, method, distance_to_similarity_list, diversity_type="class", has_similarity_matrix=False, given_similarity_matrix=None):

	if has_similarity_matrix:
		similarity_list = given_similarity_matrix

	else:
		similarity_list = get_similarity_matrix(unique_seqs_from_this_repertoires, cost, method, distance_to_similarity_list)

	zpi_list = calculate_zpi_list(similarity_list, unique_seqs_from_this_repertoires, repertoire_to_seq_prob_hash, alpha_diversity_repertoire, p)

	P_bar_dot_repertoire = repertoire_to_seq_prob_hash[alpha_diversity_repertoire]

	""" Class diversity """

	# initialize qDs
	class_alpha_diversity_results_list = defaultdict(float)
	if 1.0 in list_of_qs: class_alpha_diversity_results_list[1.] = 1. # initialize 1Ds value to 1. (defaultdict will initialize it to 0., which will cause problems with our *=, making everyting 0!)
	if float('inf') in list_of_qs: 
		class_alpha_diversity_results_list[float('inf')] = len(unique_seqs_from_this_repertoires) # initialize infDs value to the maximum possible

	for i, zpi in enumerate(zpi_list):
		for q in list_of_qs:
			if q == float('inf'): 
				class_alpha_diversity_results_list[float('inf')] = min( 1. / zpi, class_alpha_diversity_results_list[float('inf')] )
			elif q == 1.0: 
				class_alpha_diversity_results_list[q] *= 1. / ( zpi ** P_bar_dot_repertoire[i] )
			else: 
				class_alpha_diversity_results_list[q] += P_bar_dot_repertoire[i] * zpi**( q-1 )
	
	class_alpha_diversity_results_list = {q : class_alpha_diversity_results_list[q] if q in [1., float('inf')] else class_alpha_diversity_results_list[q]**( 1./(1-q) ) for q in sorted(class_alpha_diversity_results_list.keys()) }

	""" Raw diversity """
	list_of_qs_for_recon = [str(i) for i in list_of_qs]

	if unit_test:
		if "/" in recon_file_for_this_repertoire:
			recon_out_file_1 = "%s_%.2f_recon_out_1.txt" % (recon_file_for_this_repertoire.split("/")[-1].split(".")[0], cost)
			recon_out_file_2 = "%s_%.2f_recon_out_2.txt" % (recon_file_for_this_repertoire.split("/")[-1].split(".")[0], cost)
		else:
			recon_out_file_1 = "%s_%.2f_recon_out_1.txt" % (recon_file_for_this_repertoire.split(".")[0], cost)
			recon_out_file_2 = "%s_%.2f_recon_out_2.txt" % (recon_file_for_this_repertoire.split(".")[0], cost)

	else:
		if "/" in recon_file_for_this_repertoire:
			recon_out_file_1 = "%s_%.2f_recon_out_1_%s.txt" % (recon_file_for_this_repertoire.split("/")[-1].split(".")[0], cost, run_id)
			recon_out_file_2 = "%s_%.2f_recon_out_2_%s.txt" % (recon_file_for_this_repertoire.split("/")[-1].split(".")[0], cost, run_id)
		else:
			recon_out_file_1 = "%s_%.2f_recon_out_1_%s.txt" % (recon_file_for_this_repertoire.split(".")[0], cost, run_id)
			recon_out_file_2 = "%s_%.2f_recon_out_2_%s.txt" % (recon_file_for_this_repertoire.split(".")[0], cost, run_id)


	if clone_distribution_in_file: clone_distribution_option = "-c"
	else: clone_distribution_option = ""

	recon_cmd_1 = "python3 recon_v2.5.py -R %s -t 30 -l 50 -o '%s' '%s'" % (clone_distribution_option, recon_out_file_1, recon_file_for_this_repertoire)
	recon_cmd_2 = "python3 recon_v2.5.py -D -Q %s -b error_bar_parameters.txt -o %s %s" % (" ".join(list_of_qs_for_recon), recon_out_file_2, recon_out_file_1)	

	screen_out_1 = subprocess.getoutput(recon_cmd_1)
	screen_out_2 = subprocess.getoutput(recon_cmd_2)


	with open(recon_out_file_2) as f:
		for line in f:
			line = (line.strip()).split("\t")
			if '#' in line or line=='': continue
			list_of_column_name_indices = {}
			if "sample_name" in line:
				for ind, li in enumerate(line):
					list_of_column_name_indices[ind] = li
				break

	""" keys of this dict are ints, so we can arrange them """

	indices, col_names = list(zip(*sorted(list_of_column_name_indices.items())))
	est_col_names = [i for i in col_names if "est_" in i and ('+' not in i and '-' not in i)]

	with open(recon_out_file_2) as f:

		for line in f:
	
			D_number_hash = {}
			if line.startswith("#") or "sample_name" in line or line.strip() == "": continue
	
			line = (line.strip("\n")).split("\t")
	
			for ind, li in enumerate(line):
		
				try: D_number_hash[list_of_column_name_indices[ind]] = float(li)
				except: D_number_hash[list_of_column_name_indices[ind]] = li

	D_number_hash = list(D_number_hash.items())
	D_number_hash = [(i,j) for i,j in D_number_hash if "est_" in i and ('+' not in i and '-' not in i)]

	D_number_hash = [ (float(i.split("_")[1][:-1]), j) for i,j in D_number_hash ]
	raw_alpha_diversity_results_list = { i:j for i,j in D_number_hash }

	os.remove(recon_out_file_1)
	os.remove(recon_out_file_2)

	return class_alpha_diversity_results_list, raw_alpha_diversity_results_list
		

def run_beta_unit_test(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method, distance_to_similarity_list, has_similarity_matrix=False, given_similarity_matrix=None, unit_test=False):

	
	# functional_B_bar, functional_R_bar, functional_beta_bar, functional_rho_bar = generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method,  has_similarity_matrix=True, given_similarity_matrix=similarity_matrix_for_beta_diversity, unit_test=unit_test)

	functional_B_bar, functional_R_bar, functional_beta_bar, functional_rho_bar = generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method, distance_to_similarity_list, has_similarity_matrix=True, given_similarity_matrix=similarity_matrix_for_beta_diversity, unit_test=unit_test)

	print("\nTesting functional beta diversity...")
	print()

	print("\tbeta_bar")
	for reperoire_ in repertoire_names_list:
		for q in list_of_qs:
			if "%.3f" % functional_beta_bar[reperoire_][q] == '0.976': print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
			else:
				print("functional_beta_bar incorrect. Exiting...")
				exit()

	print("\n\n\trho_bar")
	for reperoire_ in repertoire_names_list:
		for q in list_of_qs:
			if "%.3f" % functional_rho_bar[reperoire_][q] == '1.025':print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
			else:
				print("functional_rho_bar incorrect. Exiting...")
				exit()
		
	print("\n\n\tB_bar")
	for q in list_of_qs:
		if "%.3f" % functional_B_bar[q] == '0.976': print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
		else:
			print("functional_B_bar incorrect. Exiting...")
			exit()

	print("\n\n\tR_bar")
	for q in list_of_qs:
		if "%.3f" % functional_R_bar[q] == '1.025': print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
		else:
			print("functional_R_bar incorrect. Exiting...")
			exit()
	print("-------------------------------------------------------")
	
	print("\n\nTesting raw beta diversity...")
	print()

	raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar = generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method, distance_to_similarity_list, diversity_type="raw")

	print("\tbeta_bar")
	for reperoire_ in repertoire_names_list:
		# for q_index, q in enumerate(list_of_qs):
		for q in list_of_qs:
			if "%.3f" % raw_beta_bar[reperoire_][q] == '2.000':print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
			else:
				print("raw_beta_bar incorrect. Exiting...")
				exit()

	print("\n\n\trho_bar")
	for reperoire_ in repertoire_names_list:
		# for q_index, q in enumerate(list_of_qs):
		for q in list_of_qs:
			if "%.3f" % raw_rho_bar[reperoire_][q] == '0.500':print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
			else:
				print("raw_rho_bar incorrect. Exiting...")
				exit()
		
	print("\n\n\tB_bar")
	for q in list_of_qs:
		if "%.3f" % raw_B_bar[q] == '2.000': print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
		else:
			print("raw_B_bar incorrect. Exiting...")
			exit()

	print("\n\n\tR_bar")
	for q in list_of_qs:
		if "%.3f" % raw_R_bar[q] == '0.500': print(("\t(%s) q=%.1f\tpass" % (reperoire_, q)))
		else:
			print("raw_R_bar incorrect. Exiting...")
			exit()
	return

def run_alpha_unit_test(input_files_list, repertoire_names_list, dir_name, list_of_qs, distance_to_similarity_list, has_similarity_matrix=False, given_similarity_matrix=None, unit_test=False):
	print("-------------------------------------------------------")
	print("\n\nTesting functional alpha diversity...")

	all_alpha_diversity_results = generate_final_alpha_diversity_output(repertoire_names_list, input_files_list, repertoire_names_list, input_files_list, list_of_qs, distance_to_similarity_list, alpha_or_beta="alpha", has_similarity_matrix=True, given_similarity_matrix=similarity_matrix_for_alpha_diversity)

	for i, (j1, j2) in list(all_alpha_diversity_results.items()):
		print(("\n\n\trepertoire\t%s" % i))
		for q in list_of_qs:
			if "%.3f" % j1[q] == '1.500': print(("\tq=%.1f\tpass" % q))
			else:
				print("functional_alpha incorrect. Exiting...")
				exit()
	print("\n\nConsult the unit test of recon for raw_alpha diversity\n")
	return


# def read_species_count_input(infile):
def read_species_count_input(input_files_list):
	d = defaultdict(float)
	non_unique_cdr3s_flag = True # initialize
	with open(input_files_list) as f:
		for line in f:
			if line.startswith('#') or line.strip() == '': continue
			try:
				cdr3, no = (line.strip()).split('\t')
			except ValueError:
				print("Line does not contain exactly 2 columns!\nYour input file is probably not in recon format (species/clone_size\tcount).\nFix this and try again. Exiting...") # Add -c option for clone size distribution as input
				exit()
			no = float(no)
			d[cdr3] += no
	return d


def read_clone_size_distribution_input(infile):
	d = defaultdict(int)
	with open(infile) as f:
		for line in f:			
			if line.startswith('#') or line.strip() == '': continue
			try: 
				clone_size, no = line.rstrip('\n').split('\t')
			except ValueError:
				print("Error: Line does not contain exactly 2 columns.\nYour input file is probably not in recon format (species/clone_size\tcount).\nFix this and try again. Exiting...") # Add -c option for clone size distribution as input
				exit()
			clone_size = int(clone_size)
			if clone_size in list(d.keys()):
				print(("Warning: detected two input lines with the same clone_size (%i).\nSumming the frequencies for both, but this is likely a sign of an error in formation of the input file" % clone_size))
			no = int(no)
			d[clone_size] += no
	return d


def collect_sequences_and_calculate_probability_terms(input_files_list, repertoire_names_list, alpha_or_beta="beta" ):
	"""
	This function gathers all unique sequences and related probability terms from input_files_list. If there > 1 file in input_files_list, this function also looks for common sequences and counts them once.

	Note: Whenever this function is called, it must be called for a single set of files.
	This means that: 
	- for beta diversity, len(repertoire_names_list) is always 2
	- for alpha diversity, len(repertoire_names_list) is always 1 (recall that even if we are calculating alpha diversity for two repertoire, each repertoire should be independent)
	"""
	seq_to_repertoire_to_count_hash = defaultdict(lambda: defaultdict(int))
	repertoire_to_total_seqs_hash = defaultdict(int)

	list_of_lists_of_seqs_from_all_repertoires = []

	for repertoire_name, input_file in zip(repertoire_names_list, input_files_list):
		if clone_distribution_in_file: 
			seq_to_count_hash = read_clone_size_distribution_input(input_file)
		else: 
			seq_to_count_hash = read_species_count_input(input_file)

		# seqs_, counts_ = zip(*seq_to_count_hash.iteritems())
		seqs_, counts_ = list(zip(*list(seq_to_count_hash.items())))

		"""Find a way to get rid of this step. This is iterating over all seqs a second time. We have already iterated over the seqs/lines when we get seq_to_count_hash. Find a way to combine both"""
		for seq_ in seqs_:
			seq_to_repertoire_to_count_hash[seq_][repertoire_name] += seq_to_count_hash[seq_]
		
		list_of_lists_of_seqs_from_all_repertoires.append(seqs_)
		repertoire_to_total_seqs_hash[repertoire_name] = sum(counts_)

	if len(list_of_lists_of_seqs_from_all_repertoires) == 2:
		unique_seqs_from_all_repertoires = []
				
		"""identify common seqs from both repertoires and put them at the very end of the the combined unqiue seqs list.
	   	We do this to make sure that we know where the common seqs are in the overall (combined) repertoire"""

		repertoire_1_seqs = sorted(list(set(list_of_lists_of_seqs_from_all_repertoires[0]))) # ensure some order
		repertoire_2_seqs = sorted(list(set(list_of_lists_of_seqs_from_all_repertoires[1]))) # ensure some order
		common_seqs = sorted(list(set(repertoire_1_seqs).intersection(set(repertoire_2_seqs))))

		repertoire_1_seqs_index = repertoire_1_seqs.index
		repertoire_2_seqs_index = repertoire_2_seqs.index

		indices_of_common_seqs_in_repertoire_1 = [ repertoire_1_seqs_index(i) for i in common_seqs ]
		indices_of_common_seqs_in_repertoire_2 = [ repertoire_2_seqs_index(i) for i in common_seqs ]

		# Remove the common seqs from individual lists, join the lists, and then append the common seqs at the very end
		unique_seqs_from_all_repertoires = [i for j, i in enumerate(repertoire_1_seqs) if j not in indices_of_common_seqs_in_repertoire_1] + [i for j, i in enumerate(repertoire_2_seqs) if j not in indices_of_common_seqs_in_repertoire_2] + common_seqs
	
	elif len(list_of_lists_of_seqs_from_all_repertoires) == 1: 
		unique_seqs_from_all_repertoires = list_of_lists_of_seqs_from_all_repertoires[0]

	else:
		print("More than two repertoires detected. Exiting...")
		exit()

	"""Get probability terms"""

	repertoire_to_seq_prob_hash = defaultdict(array) # key=repertoire; value=list of probability of (each unique seq in S) according to this repertoire. Note that there will be many zeroes because reperoires are almost disjoint
	for repertoire_name in repertoire_names_list:
		total_seqs_for_this_repertoire = float(repertoire_to_total_seqs_hash[repertoire_name])
		
		P_bar_dot_j = array([seq_to_repertoire_to_count_hash[s][repertoire_name] for s in unique_seqs_from_all_repertoires])/total_seqs_for_this_repertoire

		repertoire_to_seq_prob_hash[repertoire_name] = P_bar_dot_j

	if alpha_or_beta=="alpha":
		return unique_seqs_from_all_repertoires, repertoire_to_seq_prob_hash
	
	else: 
		w  = [1., 1.]
		p = [ dot(i, weights) for i in zip(*list(repertoire_to_seq_prob_hash.values())) ]

		return unique_seqs_from_all_repertoires, repertoire_to_seq_prob_hash, p


def get_similarity_matrix(unique_seqs_from_all_repertoires, cost, method, distance_to_similarity_list):
	"""This function returns a list of similarities between all seqs"""
	similarity_list = ([ calculate_similarity(seq_1, seq_2, cost, method, distance_to_similarity_list) for seq_2 in unique_seqs_from_all_repertoires ] for seq_1 in unique_seqs_from_all_repertoires)
	
	return similarity_list


# @jit(nopython=True)
def calculate_denominator(similarity_, repertoire_to_seq_prob_hash, repertoire_name):
	
	"""Returns the denominator of the relative uniqueness term"""

	denominator_ = dot(similarity_, repertoire_to_seq_prob_hash[repertoire_name])
	return denominator_


def calculate_zpi_list(similarity_list, unique_seqs_from_this_repertoires, repertoire_to_seq_prob_hash, repertoire_name, p, diversity_type="class"):

	"""Returns the list of zpi terms (see alpha diversity in leinster and cobbold). Used only for alpha diversity"""

	zpi_list = []

	for similarity_ in similarity_list:
		denominator_ = calculate_denominator(similarity_, repertoire_to_seq_prob_hash, repertoire_name)
		zpi_list.append(denominator_)
	
	return zpi_list


def calculate_rho_bar_term(P_bar_ij, relative_uniqueness_i, q):
	if P_bar_ij == 0.: 
		if q==1: rho_bar_term_for_i = 1. # because when p=0 x**p = 1
		else: rho_bar_term_for_i = 0. # Avoid division by 0 issues.
	else:
		if q==1: 
			rho_bar_term_for_i = relative_uniqueness_i**P_bar_ij
		else:
			rho_bar_term_for_i = P_bar_ij * ( relative_uniqueness_i**(1. - q) )

	return rho_bar_term_for_i


def calculate_rho_and_beta_bar_for_repertoire(p, repertoire_names_list, repertoire_to_seq_prob_hash, unique_seqs_from_all_repertoires, list_of_qs, weights, cost, method, distance_to_similarity_list, diversity_type="class", has_similarity_matrix=False, given_similarity_matrix=None):
	
	"""Returns a tuple of rho and beta bars for individual repertoires"""

	rho_bar_for_all_repertoires_for_all_qs, beta_bar_for_all_repertoires_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs = [ defaultdict( lambda: defaultdict(list) ) for i in range(4) ]

	for repertoire_ in repertoire_names_list:
		for q in list_of_qs:
			rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] = 0.
			beta_bar_for_all_repertoires_for_all_qs[repertoire_][q] = 0.
			weighted_rho_bar_for_all_qs[repertoire_][q] = 0.
			weighted_beta_bar_for_all_qs[repertoire_][q] = 0.

	"""
	rho_bar_for_all_repertoires_for_all_qs is the list of lists which will be returned by this function
	We will initialise this list by 0. and do the summation as we go over sequences for every repertoire, for every q
.	
	The exception here is when q=1. We cannot initialise this by 0. because here we use prod and not sum. We will initialise this by 1.

	So we will find if list_of_qs contains 1., and if so, we will initialise corresponding indices in rho_bar_for_all_repertoires_for_all_qs with 1.

	example:

	if rho_bar_for_all_repertoires_for_all_qs = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	and list_of_qs = [0., 0.5, 1.]
	then rho_bar_for_all_repertoires_for_all_qs = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]
	"""

	if 1.0 in list_of_qs:
		for repertoire_ in repertoire_names_list:
			rho_bar_for_all_repertoires_for_all_qs[repertoire_][1.0] = 1.
			beta_bar_for_all_repertoires_for_all_qs[repertoire_][1.0] = 1.

	if diversity_type=="class":
		if has_similarity_matrix: 
			print("has_similarity_matrix is TRUE")
			similarity_matrix = given_similarity_matrix
		else: similarity_matrix = get_similarity_matrix(unique_seqs_from_all_repertoires, cost, method, distance_to_similarity_list)
	
	if diversity_type == "raw":
		for ii in range(len(p)):
			numerator_ = p[ii] # similarity_matrix in raw diversity is an identity matrix
			for repertoire_ in repertoire_names_list:
				P_bar_ij = repertoire_to_seq_prob_hash[repertoire_][ii]
				denominator_ = P_bar_ij # similarity_matrix in raw diversity is an identity matrix

				relative_uniqueness_for_seq = numerator_/denominator_

				for q in list_of_qs:
					pre_sum_or_prod_rho_bar_term = calculate_rho_bar_term(P_bar_ij, relative_uniqueness_for_seq, q)
					if q==1.: 
						rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] *= pre_sum_or_prod_rho_bar_term
					else:
						rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] += pre_sum_or_prod_rho_bar_term

	elif diversity_type == "class":
		# detailed_terms=defaultdict(list)
		detailed_terms = defaultdict(lambda: defaultdict(list))
		for ii, similarity_ in enumerate(similarity_matrix):
			numerator_ = dot(similarity_, p)

			for repertoire_ in repertoire_names_list:
				P_bar_ij = repertoire_to_seq_prob_hash[repertoire_][ii]
				denominator_ = calculate_denominator(similarity_, repertoire_to_seq_prob_hash, repertoire_)

				relative_uniqueness_for_seq = numerator_/denominator_
				
				detailed_terms[repertoire_]["P_bar_ij"].append(P_bar_ij)
				detailed_terms[repertoire_]["numerator_"].append(numerator_)
				detailed_terms[repertoire_]["denominator_"].append(denominator_)
				detailed_terms[repertoire_]["relative_uniqueness_for_seq"].append(relative_uniqueness_for_seq)

				for q in list_of_qs:
					pre_sum_or_prod_rho_bar_term = calculate_rho_bar_term(P_bar_ij, relative_uniqueness_for_seq, q)
					if q==1.: 
						rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] *= pre_sum_or_prod_rho_bar_term
					else:
						rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] += pre_sum_or_prod_rho_bar_term

		# if verbose:
		# 	print("detailed_terms = ", detailed_terms)

	for repertoire_ in repertoire_names_list:
		for q in list_of_qs:
			if q != 1.:
				rho_bar_for_all_repertoires_for_all_qs[repertoire_][q] = rho_bar_for_all_repertoires_for_all_qs[repertoire_][q]**(1./(1-q))

	"""
	We double-checked that this is true: you do indeed just take the reciprocal of rho_bar, as opposed to reciprocals of the arguments of the sum. This is both explicitly stated in Table 2 in Reeve et al and calculated by Rohit to prove that Table 2 is correct

	Note that B_bar is the average of beta_bars. It is NOT 1./R_bar. R_bar is the average of rho_bars
	"""	

	for repertoire_no, repertoire_ in enumerate(repertoire_names_list):
		for q in list_of_qs:

			beta_bar_for_all_repertoires_for_all_qs[repertoire_][q] = 1./rho_bar_for_all_repertoires_for_all_qs[repertoire_][q]

			# calculate the weighted measures which are needed for B_bar and R_bar
			if q==1:
				weighted_rho_bar_for_all_qs[repertoire_][q] = (rho_bar_for_all_repertoires_for_all_qs[repertoire_][q])**weights[repertoire_no]
				weighted_beta_bar_for_all_qs[repertoire_][q] = ( beta_bar_for_all_repertoires_for_all_qs[repertoire_][q])**weights[repertoire_no]
			
			else:
				weighted_rho_bar_for_all_qs[repertoire_][q] = weights[repertoire_no]*(rho_bar_for_all_repertoires_for_all_qs[repertoire_][q])**(1.-q)
				weighted_beta_bar_for_all_qs[repertoire_][q] = weights[repertoire_no]*(beta_bar_for_all_repertoires_for_all_qs[repertoire_][q])**(1.-q)
	
	return rho_bar_for_all_repertoires_for_all_qs, beta_bar_for_all_repertoires_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs



def generate_final_alpha_diversity_output(alpha_diversity_repertoires_list, input_files_list, repertoire_names_list, recon_files_list, list_of_qs, distance_to_similarity_list, alpha_or_beta="alpha", has_similarity_matrix=False, given_similarity_matrix=None):

	all_alpha_diversity_results = {}

	""" which repertoires to calculate alpha diversity for? """

	for ii, alpha_diversity_repertoire in enumerate(alpha_diversity_repertoires_list):

		input_file_for_this_repertoire = input_files_list[  repertoire_names_list.index(alpha_diversity_repertoire) ]

		recon_file_for_this_repertoire = recon_files_list[ii]

		unique_seqs_from_this_repertoires, repertoire_to_seq_prob_hash = collect_sequences_and_calculate_probability_terms(input_file_for_this_repertoire.split(), alpha_diversity_repertoire.split(), alpha_or_beta=alpha_or_beta )

		if verbose: print(("number of unique seqs from repertoire %s:\t%i" % (alpha_diversity_repertoire, len(unique_seqs_from_this_repertoires))))

		p = np.empty(len(unique_seqs_from_this_repertoires), )

		class_alpha_diversity_results_list, raw_alpha_diversity_results_list = calculate_alpha_diversity(unique_seqs_from_this_repertoires, p, repertoire_to_seq_prob_hash, alpha_diversity_repertoire, recon_file_for_this_repertoire, list_of_qs, cost, method, distance_to_similarity_list, has_similarity_matrix=has_similarity_matrix, given_similarity_matrix=given_similarity_matrix)

		all_alpha_diversity_results[alpha_diversity_repertoire] = ( class_alpha_diversity_results_list, raw_alpha_diversity_results_list )
		
		if verbose:
			print(("alpha_diversity (%s) done:\t%s" % (alpha_diversity_repertoire, strftime("%Y-%m-%d %H:%M:%S"))))
		
		del(unique_seqs_from_this_repertoires, repertoire_to_seq_prob_hash)
	
	return all_alpha_diversity_results


def generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method, distance_to_similarity_list, diversity_type="class", has_similarity_matrix=False, given_similarity_matrix=None, unit_test=False):
	"""
	Returns all the beta diversity parameters: B_bar, R_bar, beta_bar, rho_bar, raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar
	"""
	unique_seqs_from_all_repertoires, repertoire_to_seq_prob_hash, p = collect_sequences_and_calculate_probability_terms(input_files_list, repertoire_names_list)
	
	if verbose: 
		print(("number of unique seqs from both repertoires:\t%i" % len(unique_seqs_from_all_repertoires)))

	"""Functional Diversity"""
	if verbose: print("\nfunctional beta diversity")

	rho_bar_for_all_repertoires_for_all_qs, beta_bar_for_all_repertoires_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs = calculate_rho_and_beta_bar_for_repertoire(p, repertoire_names_list, repertoire_to_seq_prob_hash, unique_seqs_from_all_repertoires, list_of_qs, weights, cost, method, distance_to_similarity_list, diversity_type=diversity_type, has_similarity_matrix=has_similarity_matrix, given_similarity_matrix=given_similarity_matrix)

	B_bar, R_bar = [ {} for i in range(2) ]
	beta_bar, rho_bar  = [ defaultdict(lambda: defaultdict(list)) for i in range(2) ]
	q_to_individual_beta_diversity_parameters_hash = defaultdict(lambda: defaultdict(list))

	for q_index, q in enumerate(list_of_qs):
		for repertoire_no, repertoire_ in enumerate(repertoire_names_list):

			q_to_individual_beta_diversity_parameters_hash[repertoire_][q].append( weighted_rho_bar_for_all_qs[repertoire_][q] )
			q_to_individual_beta_diversity_parameters_hash[repertoire_][q].append( weighted_beta_bar_for_all_qs[repertoire_][q] )

			beta_bar[repertoire_][q] = beta_bar_for_all_repertoires_for_all_qs[repertoire_][q]
			rho_bar[repertoire_][q] = rho_bar_for_all_repertoires_for_all_qs[repertoire_][q]

		repertoire_1, repertoire_2  = repertoire_names_list

		if q==1:
			R_bar_for_metarepertoire = q_to_individual_beta_diversity_parameters_hash[repertoire_1][q][0] * q_to_individual_beta_diversity_parameters_hash[repertoire_2][q][0]
			B_bar_for_metarepertoire = q_to_individual_beta_diversity_parameters_hash[repertoire_1][q][1] * q_to_individual_beta_diversity_parameters_hash[repertoire_2][q][1]
		else:	
			R_bar_for_metarepertoire = ( q_to_individual_beta_diversity_parameters_hash[repertoire_1][q][0] + q_to_individual_beta_diversity_parameters_hash[repertoire_2][q][0] )**(1./(1-q))
			B_bar_for_metarepertoire = ( q_to_individual_beta_diversity_parameters_hash[repertoire_1][q][1] + q_to_individual_beta_diversity_parameters_hash[repertoire_2][q][1] )**(1./(1-q))

		B_bar[q] = B_bar_for_metarepertoire
		R_bar[q] = R_bar_for_metarepertoire

	if verbose: print(("beta_diversity (%s) done:\t%s" % (diversity_type, strftime("%Y-%m-%d %H:%M:%S"))))
	return B_bar, R_bar, beta_bar, rho_bar




"""MAIN"""

if __name__ == '__main__':

	parser = ArgumentParser(description="""
	diversity_with_similarity: Calculate diversity with similarity for any qD
	""", formatter_class=RawDescriptionHelpFormatter)
	pa = parser.add_argument

	pa("-ar", "--alpha_diversity_repertoires", type=str, default=None, help="Alpha diversity will be calculated for both repertoires by default")
	pa("-clone_distribution_in_file", "--clone_distribution_in_file", action='store_true', help='input file is of the form clone_size tab number_of_clones of that size') # This is the same as -c option in recon
	pa("-cost", "--cost", type=float, default= 0.55, help="similarity = cost**edit_distance. Cost must be <= 1. If cost=0 then every pair of seqs is maximally different. If cost=1 then all seqs are identical") 
	pa("-dir_name", "--dir_name", type=str, default=None, help="this is the directory where all the files needed to run morty.py are stored")
	pa("-if", "--input_files", type=str, default="", help="filenames for the repertoire_names in the seq \t count format. This is usually the sampled sequences")
	pa("-ma", "--master_filename_alpha", type=str, default="alpha_diversity_master_file.txt", help="this is the master alpha output file")
	pa("-mb", "--master_filename_beta", type=str, default="beta_diversity_master_file.txt", help="this is the master beta output file")
	pa("-method", "--method", type=str, default="fast_similarity", help="method to calculate similarity. Option: fast_similarity or old_similarity")
	pa("-mo", "--mode", type=str, default="beta", help="acceptable values: alpha, beta, alpha_and_beta")
	pa("-mod", "--master_output_dir", type=str, default=None, help="this is where the master output file(s) beta_diversity_master_file.txt/alpha_diversity_master_file.txt will be written")
	pa("-qs", "--list_of_qs", type=str, default= "[0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., inf]", help="list of qs for which diversity is to be calculated")
	pa("-rn", "--repertoire_names", type=str, default="", help="names of repertoires being evaluated (>=1)")
	pa("-rf", "--recon_files", type=str, default="", help="filenames for the repertoire_names in the seq \t count format. This is usually the file with all sequences")
	pa("-sim", "--which_similarity", type=str, default= 'fast_similarity', help="choose between fast_similarity and old_similarity")
	pa("-u", "--unit_test", action="store_true", help="run unit tests")
	pa("-v", "--verbose", action="store_true", help="be verbose")
	pa("-weights", "--weights", type=str, default= "[0.5, 0.5]", help="list of weights for r1 and r2")

	args = parser.parse_args()
	globals().update(vars(args))

	#			  		#
	# *---UNIT TEST---* #
	#			  		#

	if unit_test:

		"""
		Make two files unit_test_1.txt and unit_test_2.txt

		cat unit_test_1.txt
		AAA 	1
		BBB		1
		CCC 	1

		cat unit_test_2.txt

		DDD 	1
		EEE		1
		FFF 	1
		"""

		input_files_list = ['unit_test_1.txt', 'unit_test_2.txt']
		repertoire_names_list = ['unit_test_file_1', 'unit_test_file_2']
		list_of_qs = [0., 1., 2., 3, 3.5]
		# list_of_qs = [3.]
		weights = [0.5, 0.5]
		beta_diversity=True
		
		similarity_matrix_for_beta_diversity = array([ 
			[1.,  0.5, 0.5, 0.7, 0.7, 0.7], 
			[0.5, 1.,  0.5, 0.7, 0.7, 0.7], 
			[0.5, 0.5, 1.,  0.7, 0.7, 0.7], 
			[0.7, 0.7, 0.7, 1,   0.5, 0.5],
			[0.7, 0.7, 0.7, 0.5, 1,   0.5], 
			[0.7, 0.7, 0.7, 0.5, 0.5, 1  ] ])

		similarity_matrix_for_alpha_diversity = array([
			[1.,  0.5, 0.5],
			[0.5, 1.,  0.5],
			[0.5, 0.5, 1. ] ])

		"""
		similarity_matrix_for_beta_diversity = array([ 
			[1.,  0.2, 0.2, 0.9, 0.9, 0.9], 
			[0.2, 1.,  0.2, 0.9, 0.9, 0.9], 
			[0.2, 0.2, 1.,  0.9, 0.9, 0.9], 
			[0.9, 0.9, 0.9, 1,   0.2, 0.2], 
			[0.9, 0.9, 0.9, 0.2, 1,   0.2], 
			[0.9, 0.9, 0.9, 0.2, 0.2, 1  ] ])

		similarity_matrix_for_beta_diversity = array([ 
			[1.,  0., 0., 1., 1., 1.], 
			[0., 1.,  0., 1., 1., 1.], 
			[0., 0., 1.,  1., 1., 1.], 
			[1., 1., 1., 1,   0., 0.], 
			[1., 1., 1., 0., 1,   0.], 
			[1., 1., 1., 0., 0., 1  ] ])
		"""

		if not dir_name.endswith("/"): dir_name = dir_name+"/"

		distance_to_similarity_list = get_distance_to_similarity_list(cost, method)

		run_beta_unit_test(input_files_list, repertoire_names_list, dir_name, list_of_qs, cost, method, distance_to_similarity_list, has_similarity_matrix=True, given_similarity_matrix=similarity_matrix_for_beta_diversity, unit_test=unit_test)

		run_alpha_unit_test(input_files_list, repertoire_names_list, dir_name, list_of_qs, distance_to_similarity_list, has_similarity_matrix=True, given_similarity_matrix=similarity_matrix_for_alpha_diversity, unit_test=unit_test)

		exit()

	date_ = strftime("%Y-%m-%d %H:%M:%S")
	print(); print((strftime("%Y-%m-%d %H:%M:%S"))); print()
	start_time = time()

	alpha_diversity=False
	beta_diversity=False

	if mode=="alpha": alpha_diversity=True
	elif mode=="beta": beta_diversity=True
	elif mode=="alpha_and_beta":
		alpha_diversity=True
		beta_diversity=True
	else:
		print("Specify mode. Exiting...")
		exit()

	# Set output files
	if not master_output_dir:

		master_output_dirs = [expanduser("~") + "/Dropbox (ArnaoutLab)/Harry/repo/Python scripts/", "/n/data2/bidmc/path/arnaout/", expanduser("~") + "/Desktop/Morty"]
		master_output_dir_which_exist = [i for i in master_output_dirs if os.path.isdir(i)]
		# master_output_dir = master_output_dir_which_exist[0]
		if len(master_output_dir_which_exist) == 0:
			print("Provide path to master output files. Exiting...")
			exit()

	master_filename_beta = "%s/%s" % (master_output_dir, master_filename_beta)
	master_filename_alpha = "%s/%s" % (master_output_dir, master_filename_alpha)

	recon_files_list = recon_files.split(',')
	repertoire_names_list = repertoire_names.split(',')
	input_files_list = input_files.split(',')

	set_cython_seed(randint(0, 2147483647)) # This is cython RAND_MAX

	code_version = argv[0]
	datestamp = subprocess.getoutput("date +%m%d%Y")

	if not dir_name.endswith("/"): dir_name = dir_name+"/"

	run_id = ''.join(choice(ascii_letters + digits) for i in range(5))
	if verbose: print(("run_id = %s" % run_id))

	results_dir = "%sBeta_diversity_results" % dir_name
	if not os.path.isdir(results_dir): os.mkdir(results_dir)

	# Process list of qs
	list_of_qs = list_of_qs.replace("inf", "float('inf')")
	list_of_qs = eval(list_of_qs)

	weights = eval(weights)

	recon_file_to_repertoire_name = list(zip(repertoire_names_list, input_files_list))

	if not alpha_diversity_repertoires: alpha_diversity_repertoires_list = repertoire_names_list
	else: alpha_diversity_repertoires_list = alpha_diversity_repertoires.split()

	"""get column names for d numbers"""
	col_names_for_qs = []
	for q in [str(i) for i in list_of_qs]:
		col_name = "%sDs %sD" % (q, q)
		col_names_for_qs.append(col_name)

	col_names_for_qs = list(zip(*[i.split() for i  in col_names_for_qs]))

	screen_out = ''; outstr_beta_measures=''; outstr_alpha_measures=''

	if not os.path.isfile(master_filename_beta):
		outstr_beta_measures += "run_id\tmeasure\trep_1\trep_2\tq_hash\tdate\tcommand\n"
	if not os.path.isfile(master_filename_alpha):
		outstr_alpha_measures += "run_id\tmeasure\trep\tq_hash\tdate\tcommand\n"
		
	if verbose:
		for i in recon_file_to_repertoire_name:
			i1, i2 = i
			screen_out += "# Input for %s:\t%s\n" % (i1, i2)

	if alpha_diversity: screen_out += "\n# ALPHA DIVERSITY\n# ---------------\n# status\t%s\n# repertoires\t%s\n\n\n" % ( alpha_diversity, ", ".join(alpha_diversity_repertoires_list) )
	if beta_diversity: screen_out += "# BETA DIVERSITY\n# ---------------\n# status\t%s\n# repertoires\t%s\n\n" % ( beta_diversity, ", ".join(repertoire_names_list) )

	if verbose: print(screen_out)

	if alpha_diversity:
		if verbose: print("Calculating alpha diversity")

		all_alpha_diversity_results = generate_final_alpha_diversity_output(alpha_diversity_repertoires_list, input_files_list, repertoire_names_list, recon_files_list, list_of_qs, alpha_or_beta="alpha")


		outstr_alpha = ''
		for repertoire_ in alpha_diversity_repertoires_list:

			functional_alpha_diversity_results, raw_alpha_diversity_results = all_alpha_diversity_results[repertoire_]
				
			functional_alpha_diversity_results = { str(q)+"Ds" : "%.3e"%functional_alpha_diversity_results[q] for q in list(functional_alpha_diversity_results.keys()) }
			raw_alpha_diversity_results = { str(q)+"Ds" : "%.3e"%raw_alpha_diversity_results[q] for q in list(raw_alpha_diversity_results.keys()) }

			outstr_alpha += "%s\talpha\t%s\t" % (run_id, repertoire_)
			class_vals = "(%s" % functional_alpha_diversity_results
			outstr_alpha += "%s, " % class_vals
			raw_vals = "%s)\t" % raw_alpha_diversity_results
			outstr_alpha += "%s" % raw_vals
			outstr_alpha += "%s\tpython %s\n" % (date_, " ".join(argv))

		outstr_alpha_measures += outstr_alpha

		if verbose: print(("\n\n%s" % outstr_alpha_measures))

		if os.path.isfile(master_filename_alpha): action_ = 'a'
		else: action_ = 'w'
		with open(master_filename_alpha, action_) as f: f.write(outstr_alpha_measures)
		
	if beta_diversity:
		if verbose: 
			print("\nCalculating beta diversity")
		
		if float('inf') in list_of_qs:
			list_of_qs_beta = [ q for q in list_of_qs if q != float("inf") ]
		elif float('-inf') in list_of_qs:
			list_of_qs_beta = [ q for q in list_of_qs if q != float("-inf") ]
		else: 
			list_of_qs_beta = list_of_qs
		
		distance_to_similarity_list = get_distance_to_similarity_list(cost, method)
		
		functional_B_bar, functional_R_bar, functional_beta_bar, functional_rho_bar = generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs_beta, cost, method, distance_to_similarity_list)

		functional_B_bar = { str(q)+"Ds" : "%.3e"%functional_B_bar[q] for q in list(functional_B_bar.keys()) }
		functional_R_bar = { str(q)+"Ds" : "%.3e"%functional_R_bar[q] for q in list(functional_R_bar.keys()) }

		raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar = generate_final_beta_diversity_output(input_files_list, repertoire_names_list, dir_name, list_of_qs_beta, cost, method, distance_to_similarity_list, diversity_type="raw")

		raw_B_bar = {str(q)+"D" : "%.3e"%raw_B_bar[q] for q in list(raw_B_bar.keys())}
		raw_R_bar = {str(q)+"D" : "%.3e"%raw_R_bar[q] for q in list(raw_R_bar.keys())}


		""" first do B_bar and R_bar """
		outstr_B = ''; outstr_R = ''

		repertoire_1, repertoire_2 = repertoire_names_list

		if float('inf') in list_of_qs:
			raw_B_bar['infD'] = float('nan')
			raw_R_bar['infD'] = float('nan')
			functional_B_bar['infDs'] = float('nan')
			functional_R_bar['infDs'] = float('nan')

		if float('-inf') in list_of_qs:
			raw_B_bar['-infD'] = float('nan')
			raw_R_bar['-infD'] = float('nan')
			functional_B_bar['-infDs'] = float('nan')
			functional_R_bar['-infDs'] = float('nan')

		functional_raw_B_bar = dict(list(functional_B_bar.items()) + list(raw_B_bar.items()))
		functional_raw_R_bar = dict(list(functional_R_bar.items()) + list(raw_R_bar.items()))

		outstr_B += "%s\tB_bar\t%s\t%s\t" % (run_id, repertoire_1, repertoire_2)
		outstr_B += "%s\t" % str(functional_raw_B_bar)
		outstr_B += "%s\tpython %s\n" % (date_, " ".join(argv))
		outstr_beta_measures += outstr_B

		outstr_R += "%s\tR_bar\t%s\t%s\t" % (run_id, repertoire_1, repertoire_2)
		outstr_R += "%s\t" % str(functional_raw_R_bar)
		outstr_R += "%s\tpython %s\n" % (date_, " ".join(argv))
		outstr_beta_measures += outstr_R
		
		""" Now do beta_bar and rho_bar """
		for repertoire_ in repertoire_names_list:
			other_repertoire = [i for i in repertoire_names_list if i != repertoire_][0]
			functional_beta_bar[repertoire_] = { str(q)+"Ds" : "%.3e"%functional_beta_bar[repertoire_][q] for q in list(functional_beta_bar[repertoire_].keys()) }
			raw_beta_bar[repertoire_] = { str(q)+"D" : "%.3e"%raw_beta_bar[repertoire_][q] for q in list(raw_beta_bar[repertoire_].keys()) }

			outstr_beta = ''
			if float('inf') in list_of_qs:
				functional_beta_bar[repertoire_]['infDs'] = float('nan')
				raw_beta_bar[repertoire_]['infD'] = float('nan')
			#
			if float('-inf') in list_of_qs:
				functional_beta_bar[repertoire_]['-infDs'] = float('nan')
				raw_beta_bar[repertoire_]['-infD']  = float('nan')

			functional_raw_beta_bar = {}
			for k_ in list(functional_beta_bar.keys()):
				functional_raw_beta_bar[k_] = dict( list(functional_beta_bar[k_].items()) + list(raw_beta_bar[k_].items()) )

			outstr_beta += "%s\tbeta_bar\t%s\t%s\t" % (run_id, repertoire_, other_repertoire)
			outstr_beta += "%s\t" % dict(functional_raw_beta_bar[repertoire_])
			outstr_beta += "%s\tpython %s\n" % (date_, " ".join(argv))
			outstr_beta_measures += outstr_beta

		for repertoire_ in repertoire_names_list:
			
			other_repertoire = [i for i in repertoire_names_list if i != repertoire_][0]
			
			functional_rho_bar[repertoire_] = { str(q)+"Ds" : "%.3e"%functional_rho_bar[repertoire_][q] for q in list(functional_rho_bar[repertoire_].keys()) }
			raw_rho_bar[repertoire_] = { str(q)+"D" : "%.3e"%raw_rho_bar[repertoire_][q] for q in list(raw_rho_bar[repertoire_].keys()) }

			outstr_rho = ''
			if float('inf') in list_of_qs:
				functional_rho_bar[repertoire_]['infDs'] = float('nan')
				raw_rho_bar[repertoire_]['infD']  = float('nan')
			#
			if float('-inf') in list_of_qs:
				functional_rho_bar[repertoire_]['-infDs'] = float('nan')
				raw_rho_bar[repertoire_]['-infD']  = float('nan')

			functional_raw_rho_bar = {}
			for k_ in list(functional_rho_bar.keys()):
				functional_raw_rho_bar[k_] = dict( list(functional_rho_bar[k_].items()) + list(raw_rho_bar[k_].items()) )

			outstr_rho += "%s\trho_bar\t%s\t%s\t" % (run_id, repertoire_, other_repertoire)
			outstr_rho += "%s\t" % dict(functional_raw_rho_bar[repertoire_])
			outstr_rho += "%s\tpython %s\n" % (date_, " ".join(argv))
			outstr_beta_measures += outstr_rho

		if verbose: print(("\n\n%s" % outstr_beta_measures))

		if os.path.isfile(master_filename_beta): action_ = 'a'
		else: action_ = 'w'
		with open(master_filename_beta, action_) as f: f.write(outstr_beta_measures)

	print(); print((strftime("%Y-%m-%d %H:%M:%S")))
	if verbose: print(("# time_taken\t%.2fs" % (time() - start_time)))