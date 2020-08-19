#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser, RawDescriptionHelpFormatter, RawTextHelpFormatter
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
from simlib import set_cython_seed
from string import *
from sys import argv, exit
from time import time, sleep, strftime
import numpy as np
import os
import subprocess
import importlib
import sys
import warnings


warnings.filterwarnings("ignore", category=RuntimeWarning) 
path_to_morty = "/".join(subprocess.getoutput("which morty.py").split("/")[:-1])

""" FUNCTIONS """

def eval_right_side(txt):
	return literal_eval(txt.split("=")[1].str.strip())


def prod(iterable):
	"""The equivalent of sum() for computing the product of the terms
	Arguments:
		iterable: takes an iterable [such as a list] and returns the product of the values"""
	return reduce(mul, iterable, 1)


def read_species_count_input(input_files_list):
	d = defaultdict(float)
	non_unique_cdr3s_flag = True # initialize
	with open(input_files_list) as f:
		for line in f:
			if line.startswith('#') or line.strip() == '': continue
			if 'cdr3_aa' in line: continue
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
			if 'cdr3_aa' in line: continue
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


def parse_input_files(input_files_list, community_names_list):
	"""
	This function gathers all unique sequences and related probability terms from input_files_list. If there > 1 file in input_files_list, this function also looks for common sequences and counts them once.

	Note: Whenever this function is called, it must be called for a single set of files.
	This means that: 
	- for beta diversity, len(community_names_list) is always 2
	- for alpha diversity, len(community_names_list) is always 1 (recall that even if we are calculating alpha diversity for two repertoire, each repertoire should be independent)
	"""
	
	species_to_community_to_count_hash = defaultdict(lambda: defaultdict(int))
	community_to_total_species_hash = defaultdict(int)

	list_of_lists_of_seqs_from_all_repertoires = []

	for community_name, input_file in zip(community_names_list, input_files_list):
		if clone_distribution_in_file: 
			seq_to_count_hash = read_clone_size_distribution_input(input_file)
		else: 
			seq_to_count_hash = read_species_count_input(input_file)

		seqs_, counts_ = list(zip(*list(seq_to_count_hash.items())))

		"""Find a way to get rid of this step. This is iterating over all seqs a second time. We have already iterated over the seqs/lines when we get seq_to_count_hash. Find a way to combine both"""
		for seq_ in seqs_:
			species_to_community_to_count_hash[seq_][community_name] += seq_to_count_hash[seq_]
		
		list_of_lists_of_seqs_from_all_repertoires.append(seqs_)
		community_to_total_species_hash[community_name] = sum(counts_)

	return list_of_lists_of_seqs_from_all_repertoires, community_to_total_species_hash, species_to_community_to_count_hash


def get_unique_sequences_from_input_files(list_of_lists_of_seqs_from_all_repertoires):

	if len(list_of_lists_of_seqs_from_all_repertoires) == 2:
		unique_species_from_all_communities = []
		
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
		unique_species_from_all_communities = [i for j, i in enumerate(repertoire_1_seqs) if j not in indices_of_common_seqs_in_repertoire_1] + [i for j, i in enumerate(repertoire_2_seqs) if j not in indices_of_common_seqs_in_repertoire_2] + common_seqs
	
	elif len(list_of_lists_of_seqs_from_all_repertoires) == 1: 
		unique_species_from_all_communities = list_of_lists_of_seqs_from_all_repertoires[0]

	else:
		print("More than two communities detected.\nExiting...")
		exit()
	
	return unique_species_from_all_communities


def get_seqs_from_unique_seqs_user(unique_species_user):
	unique_seqs_user_list=[]
	with open(unique_species_user) as f:
		for line in f:
			line = str.strip(line, "\n")
			unique_seqs_user_list.append(line)
	return unique_seqs_user_list


def calculate_probability_terms(unique_species_from_all_communities, community_names_list, community_to_total_species_hash, species_to_community_to_count_hash, alpha_or_beta="beta") :
	"""Get probability terms"""
	community_to_species_prob_hash = defaultdict(array) # key=repertoire; value=list of probability of (each unique seq in S) according to this repertoire. Note that there will be many zeroes because reperoires are almost disjoint
	for community_name in community_names_list:
		total_species_for_this_community = float(community_to_total_species_hash[community_name])
		
		P_bar_dot_j = array([species_to_community_to_count_hash[s][community_name] for s in unique_species_from_all_communities]) / total_species_for_this_community

		community_to_species_prob_hash[community_name] = P_bar_dot_j

	if alpha_or_beta=="alpha":
		return community_to_species_prob_hash
	
	else: 
		w  = [1., 1.]
		p = [ dot(i, weights) for i in zip(*list(community_to_species_prob_hash.values())) ]

		return community_to_species_prob_hash, p


def calculate_alpha_diversity(unique_species, p, community_to_species_prob_hash, community_name, recon_file, list_of_qs, diversity_type="class", user_similarity_matrix=None, user_function=None):

	community_name = community_name[0] # recall that community_name in the main function is a lisr

	if user_similarity_matrix is not None: 
		similarity_matrix = user_similarity_matrix
	
	elif user_function:
		similarity_matrix = user_function(unique_species)

	zpi_list = calculate_zpi_list(similarity_matrix, unique_species, community_to_species_prob_hash, community_name, p)
	
	P_bar_dot_community = community_to_species_prob_hash[community_name]

	""" Class diversity """
	# initialize qDs
	class_alpha_diversity_results_list = defaultdict(float)
	if 1.0 in list_of_qs: class_alpha_diversity_results_list[1.] = 1. # initialize 1Ds value to 1. (defaultdict will initialize it to 0., which will cause problems with our *=, making everyting 0!)
	if float('inf') in list_of_qs: 
		class_alpha_diversity_results_list[float('inf')] = len(unique_species) # initialize infDs value to the maximum possible

	for i, zpi in enumerate(zpi_list):
		for q in list_of_qs:
			if q == float('inf'): 
				class_alpha_diversity_results_list[float('inf')] = min( 1. / zpi, class_alpha_diversity_results_list[float('inf')] )
			elif q == 1.0: 
				class_alpha_diversity_results_list[q] *= 1. / ( zpi ** P_bar_dot_community[i] )
			else: 
				class_alpha_diversity_results_list[q] += P_bar_dot_community[i] * zpi**( q-1 )
	
	class_alpha_diversity_results_list = {q : class_alpha_diversity_results_list[q] if q in [1., float('inf')] else class_alpha_diversity_results_list[q]**( 1./(1-q) ) for q in sorted(class_alpha_diversity_results_list.keys()) }

	""" Raw diversity """
	list_of_qs_for_recon = [str(i) for i in list_of_qs]

	if unit_test:
		if "/" in recon_file:
			recon_out_file_1 = "%s_recon_out_1.txt" % recon_file.split("/")[-1].split(".")[0]
			recon_out_file_2 = "%s_recon_out_2.txt" % recon_file.split("/")[-1].split(".")[0]
		else:
			recon_out_file_1 = "%s_recon_out_1.txt" % recon_file.split(".")[0]
			recon_out_file_2 = "%s_recon_out_2.txt" % recon_file.split(".")[0]

	else:
		if "/" in recon_file:
			recon_out_file_1 = "%s_recon_out_1_%s.txt" % (recon_file.split("/")[-1].split(".")[0], run_id)
			recon_out_file_2 = "%s_recon_out_2_%s.txt" % (recon_file.split("/")[-1].split(".")[0], run_id)
		else:
			recon_out_file_1 = "%s_recon_out_1_%s.txt" % (recon_file.split(".")[0], run_id)
			recon_out_file_2 = "%s_recon_out_2_%s.txt" % (recon_file.split(".")[0], run_id)


	if clone_distribution_in_file: clone_distribution_option = "-c"
	else: clone_distribution_option = ""

	recon_cmd_1 = "recon_v3.0.py -R %s -t 30 -l 50 -o '%s' '%s'" % (clone_distribution_option, recon_out_file_1, recon_file)

	which_recon = subprocess.getoutput("which recon_v3.0.py")
	recon_path, _ = os.path.split(which_recon)

	recon_cmd_2 = "recon_v3.0.py -D -Q %s -b '%s'/error_bar_parameters.txt -o %s %s" % (" ".join(list_of_qs_for_recon), recon_path, recon_out_file_2, recon_out_file_1)

	screen_out_1 = subprocess.getoutput(recon_cmd_1)

	"""If recon is not found then print the error message and exit gracefully"""
	if "command not found" in screen_out_1:
		error_message = ''
		error_message += "\nRecon is either not installed or not in user's PATH. To fix this:\n"
		error_message += "i\tDownload Recon from here: https://github.com/ArnaoutLab/Recon\n"
		error_message += "ii\tAdd it your user PATH"
		print(error_message)
		exit()
	
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

	""" Keys of this dict are ints, so we can arrange them """
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
		

def calculate_denominator(similarity_, community_to_species_prob_hash, community_name):
	
	"""Returns the denominator of the relative uniqueness term"""
	denominator_ = dot(similarity_, community_to_species_prob_hash[community_name])
	return denominator_


def calculate_zpi_list(similarity_matrix, unique_species, community_to_species_prob_hash, community_name, p, diversity_type="class"):

	"""Returns the list of zpi terms (see alpha diversity in leinster and cobbold). Used only for alpha diversity"""

	zpi_list = []

	for similarity_ in similarity_matrix:
		denominator_ = calculate_denominator(similarity_, community_to_species_prob_hash, community_name)
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


def calculate_rho_and_beta_bar_for_repertoire(p, community_names_list, community_to_species_prob_hash, unique_species_from_all_communities, list_of_qs, weights, diversity_type="class", user_similarity_matrix=None, user_function=None):
	
	"""Returns a tuple of rho and beta bars for individual repertoires"""

	rho_bar_for_all_communities_for_all_qs, beta_bar_for_all_communities_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs = [ defaultdict( lambda: defaultdict(list) ) for i in range(4) ]

	for community_ in community_names_list:
		for q in list_of_qs:
			rho_bar_for_all_communities_for_all_qs[community_][q] = 0.
			beta_bar_for_all_communities_for_all_qs[community_][q] = 0.
			weighted_rho_bar_for_all_qs[community_][q] = 0.
			weighted_beta_bar_for_all_qs[community_][q] = 0.

	"""rho_bar_for_all_communities_for_all_qs is the list of lists which will be returned by this function
	We will initialise this list by 0. and do the summation as we go over sequences for every repertoire, for every q
.	
	The exception here is when q=1. We cannot initialise this by 0. because here we use prod and not sum. We will initialise this by 1.

	So we will find if list_of_qs contains 1., and if so, we will initialise corresponding indices in rho_bar_for_all_communities_for_all_qs with 1.

	Example:

	if rho_bar_for_all_communities_for_all_qs = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	and list_of_qs = [0., 0.5, 1.]
	then rho_bar_for_all_communities_for_all_qs = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]"""

	if 1.0 in list_of_qs:
		for community_ in community_names_list:
			rho_bar_for_all_communities_for_all_qs[community_][1.0] = 1.
			beta_bar_for_all_communities_for_all_qs[community_][1.0] = 1.

	if diversity_type=="class":
		"""Calculate or use custom similarity matrix"""

		if user_similarity_matrix is not None:
			similarity_matrix = user_similarity_matrix
		
		elif user_function:
			similarity_matrix = user_function(unique_species_from_all_communities)

	if diversity_type == "raw":
		for ii in range(len(p)):
			numerator_ = p[ii] # similarity_matrix in raw diversity is an identity matrix
			for community_ in community_names_list:
				P_bar_ij = community_to_species_prob_hash[community_][ii]
				denominator_ = P_bar_ij # similarity_matrix in raw diversity is an identity matrix

				relative_uniqueness_for_seq = numerator_/denominator_

				for q in list_of_qs:
					pre_sum_or_prod_rho_bar_term = calculate_rho_bar_term(P_bar_ij, relative_uniqueness_for_seq, q)
					if q==1.: 
						rho_bar_for_all_communities_for_all_qs[community_][q] *= pre_sum_or_prod_rho_bar_term
					else:
						rho_bar_for_all_communities_for_all_qs[community_][q] += pre_sum_or_prod_rho_bar_term

	elif diversity_type == "class":
		# print("similarity_matrix =", [i for i in similarity_matrix])
		detailed_terms = defaultdict(lambda: defaultdict(list))
		for ii, similarity_ in enumerate(similarity_matrix):
			numerator_ = dot(similarity_, p)

			for community_ in community_names_list:
				P_bar_ij = community_to_species_prob_hash[community_][ii]
				denominator_ = calculate_denominator(similarity_, community_to_species_prob_hash, community_)

				relative_uniqueness_for_seq = numerator_/denominator_
				
				detailed_terms[community_]["P_bar_ij"].append(P_bar_ij)
				detailed_terms[community_]["numerator_"].append(numerator_)
				detailed_terms[community_]["denominator_"].append(denominator_)
				detailed_terms[community_]["relative_uniqueness_for_seq"].append(relative_uniqueness_for_seq)

				for q in list_of_qs:
					pre_sum_or_prod_rho_bar_term = calculate_rho_bar_term(P_bar_ij, relative_uniqueness_for_seq, q)
					if q==1.: 
						rho_bar_for_all_communities_for_all_qs[community_][q] *= pre_sum_or_prod_rho_bar_term
					else:
						rho_bar_for_all_communities_for_all_qs[community_][q] += pre_sum_or_prod_rho_bar_term

	for community_ in community_names_list:
		for q in list_of_qs:
			if q != 1.:
				rho_bar_for_all_communities_for_all_qs[community_][q] = rho_bar_for_all_communities_for_all_qs[community_][q]**(1./(1-q))

	"""We double-checked that this is true: you do indeed just take the reciprocal of rho_bar, as opposed to reciprocals of the arguments of the sum. This is both explicitly stated in Table 2 in Reeve et al and calculated by Rohit to prove that Table 2 is correct

	Note that B_bar is the average of beta_bars. It is NOT 1./R_bar. R_bar is the average of rho_bars"""	

	for community_no, community_ in enumerate(community_names_list):
		for q in list_of_qs:

			beta_bar_for_all_communities_for_all_qs[community_][q] = 1./rho_bar_for_all_communities_for_all_qs[community_][q]

			# calculate the weighted measures which are needed for B_bar and R_bar
			if q==1:
				weighted_rho_bar_for_all_qs[community_][q] = (rho_bar_for_all_communities_for_all_qs[community_][q])**weights[community_no]
				weighted_beta_bar_for_all_qs[community_][q] = ( beta_bar_for_all_communities_for_all_qs[community_][q])**weights[community_no]
			
			else:
				weighted_rho_bar_for_all_qs[community_][q] = weights[community_no]*(rho_bar_for_all_communities_for_all_qs[community_][q])**(1.-q)
				weighted_beta_bar_for_all_qs[community_][q] = weights[community_no]*(beta_bar_for_all_communities_for_all_qs[community_][q])**(1.-q)
	
	return rho_bar_for_all_communities_for_all_qs, beta_bar_for_all_communities_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs


def generate_final_alpha_diversity_output(input_file, community_name, recon_file, list_of_qs, user_similarity_matrix=None, unique_species_user=None, user_function=None):

	recon_file=recon_file[0]
	all_alpha_diversity_results = {}

	list_of_lists_of_species, community_to_total_species_hash, species_to_community_to_count_hash = parse_input_files(input_file, community_name)

	if user_similarity_matrix is not None:
		unique_species = get_seqs_from_unique_seqs_user(unique_species_user)
	else: 
		unique_species = get_unique_sequences_from_input_files(list_of_lists_of_species)

	community_to_species_prob_hash = calculate_probability_terms(unique_species, community_name, community_to_total_species_hash, species_to_community_to_count_hash, alpha_or_beta="alpha")

	"""Assert that the sequences in unique_species and in the input communities are the same."""
	if user_similarity_matrix is not None:
		try:
			species_collected_from_input_files = [ j for i in list_of_lists_of_species for j in i ]
			assert sorted(unique_species) == sorted(list(set(species_collected_from_input_files)))
			
		except AssertionError:
			print("Species list in the %s and input files (%s) should be identical. See manual for details and examples.\nExiting..." % (unique_species_user, community_name))
			exit()

	if verbose: print("# Number of unique species from community %s:\t%i" % (community_name[0], len(unique_species)))

	p = np.empty(len(unique_species), ) # dummy p

	class_alpha_diversity_results_list, raw_alpha_diversity_results_list = calculate_alpha_diversity(unique_species, p, community_to_species_prob_hash, community_name, recon_file, list_of_qs, user_similarity_matrix=user_similarity_matrix, user_function=user_function)

	all_alpha_diversity_results[community_name[0]] = (class_alpha_diversity_results_list, raw_alpha_diversity_results_list)
	
	# if verbose: print(("alpha_diversity (%s) done:\t%s" % (community_name[0], strftime("%Y-%m-%d %H:%M:%S"))))
	del(unique_species, community_to_species_prob_hash)
	
	return all_alpha_diversity_results


def generate_final_beta_diversity_output(input_files_list, community_names_list, list_of_qs, diversity_type="class", user_similarity_matrix=None, user_function=None, unique_species_user=None):
	
	"""Returns all the beta diversity parameters: B_bar, R_bar, beta_bar, rho_bar, raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar"""
	
	list_of_lists_of_species, community_to_total_species_hash, species_to_community_to_count_hash = parse_input_files(input_files_list, community_names_list)

	if user_similarity_matrix is not None: 
		unique_species_user = get_seqs_from_unique_seqs_user(unique_species_user)
		unique_species_from_all_communities = unique_species_user

	else: 
		unique_species_from_all_communities = get_unique_sequences_from_input_files(list_of_lists_of_species)
	
	community_to_species_prob_hash, p = calculate_probability_terms(unique_species_from_all_communities, community_names_list, community_to_total_species_hash, species_to_community_to_count_hash)
	
	if diversity_type=="class":
		if verbose: print(("Number of unique species from both communities:\t%i" % len(unique_species_from_all_communities)))

	"""Functional Diversity"""
	rho_bar_for_all_communities_for_all_qs, beta_bar_for_all_communities_for_all_qs, weighted_rho_bar_for_all_qs, weighted_beta_bar_for_all_qs = calculate_rho_and_beta_bar_for_repertoire(p, community_names_list, community_to_species_prob_hash, unique_species_from_all_communities, list_of_qs, weights,diversity_type=diversity_type, user_similarity_matrix=user_similarity_matrix, user_function=user_function)

	B_bar, R_bar = [ {} for i in range(2) ]
	beta_bar, rho_bar  = [ defaultdict(lambda: defaultdict(list)) for i in range(2) ]
	q_to_individual_beta_diversity_parameters_hash = defaultdict(lambda: defaultdict(list))

	for q_index, q in enumerate(list_of_qs):
		for community_no, community_ in enumerate(community_names_list):

			q_to_individual_beta_diversity_parameters_hash[community_][q].append( weighted_rho_bar_for_all_qs[community_][q] )
			q_to_individual_beta_diversity_parameters_hash[community_][q].append( weighted_beta_bar_for_all_qs[community_][q] )

			beta_bar[community_][q] = beta_bar_for_all_communities_for_all_qs[community_][q]
			rho_bar[community_][q] = rho_bar_for_all_communities_for_all_qs[community_][q]

		community_1, community_2  = community_names_list

		if q==1:
			R_bar_for_metacommunity = q_to_individual_beta_diversity_parameters_hash[community_1][q][0] * q_to_individual_beta_diversity_parameters_hash[community_1][q][0]
			B_bar_for_metacommunity = q_to_individual_beta_diversity_parameters_hash[community_1][q][1] * q_to_individual_beta_diversity_parameters_hash[community_2][q][1]
		else:	
			R_bar_for_metacommunity = ( q_to_individual_beta_diversity_parameters_hash[community_1][q][0] + q_to_individual_beta_diversity_parameters_hash[community_2][q][0] )**(1./(1-q))
			B_bar_for_metacommunity = ( q_to_individual_beta_diversity_parameters_hash[community_1][q][1] + q_to_individual_beta_diversity_parameters_hash[community_2][q][1] )**(1./(1-q))

		B_bar[q] = B_bar_for_metacommunity
		R_bar[q] = R_bar_for_metacommunity

	if verbose: print(("beta_diversity (%s) done:\t%s" % (diversity_type, strftime("%Y-%m-%d %H:%M:%S"))))
	return B_bar, R_bar, beta_bar, rho_bar


def run_beta_unit_test(input_files_list, community_names_list, list_of_qs, similarity_matrix_for_beta_diversity, unique_species_user, unit_test=False):

	functional_B_bar, functional_R_bar, functional_beta_bar, functional_rho_bar = generate_final_beta_diversity_output(input_files_list, community_names_list, list_of_qs, user_similarity_matrix=similarity_matrix_for_beta_diversity, unique_species_user=unique_species_user)

	print("\nTesting functional beta diversity...\n")

	print("\tbeta_bar")
	for community_ in community_names_list:
		for q in list_of_qs:
			if "%.3f" % functional_beta_bar[community_][q] == '0.976': print(("\t(%s) q=%.1f\tpass" % (community_, q)))
			else:
				print("functional_beta_bar incorrect.\nExiting...")
				exit()

	print("\n\n\trho_bar")
	for community_ in community_names_list:
		for q in list_of_qs:
			if "%.3f" % functional_rho_bar[community_][q] == '1.025':print(("\t(%s) q=%.1f\tpass" % (community_, q)))
			else:
				print("functional_rho_bar incorrect. Exiting...")
				exit()
		
	print("\n\n\tB_bar")
	for q in list_of_qs:
		if "%.3f" % functional_B_bar[q] == '0.976': print(("\t(%s) q=%.1f\tpass" % (community_, q)))
		else:
			print("functional_B_bar incorrect. Exiting...")
			exit()

	print("\n\n\tR_bar")
	for q in list_of_qs:
		if "%.3f" % functional_R_bar[q] == '1.025': print(("\t(%s) q=%.1f\tpass" % (community_, q)))
		else:
			print("functional_R_bar incorrect. Exiting...")
			exit()
	print("-------------------------------------------------------")
	
	print("\n\nTesting raw beta diversity...")
	print("")

	raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar = generate_final_beta_diversity_output(input_files_list, community_names_list, list_of_qs, diversity_type="raw")

	print("\tbeta_bar")
	for community_ in community_names_list:
		# for q_index, q in enumerate(list_of_qs):
		for q in list_of_qs:
			if "%.3f" % raw_beta_bar[community_][q] == '2.000':print(("\t(%s) q=%.1f\tpass" % (community_, q)))
			else:
				print("raw_beta_bar incorrect. Exiting...")
				exit()

	print("\n\n\trho_bar")
	for community_ in community_names_list:
		for q in list_of_qs:
			if "%.3f" % raw_rho_bar[community_][q] == '0.500':print(("\t(%s) q=%.1f\tpass" % (community_, q)))
			else:
				print("raw_rho_bar incorrect. Exiting...")
				exit()
		
	print("\n\n\tB_bar")
	for q in list_of_qs:
		if "%.3f" % raw_B_bar[q] == '2.000': print(("\t(%s) q=%.1f\tpass" % (community_, q)))
		else:
			print("raw_B_bar incorrect. Exiting...")
			exit()

	print("\n\n\tR_bar")
	for q in list_of_qs:
		if "%.3f" % raw_R_bar[q] == '0.500': print(("\t(%s) q=%.1f\tpass" % (community_, q)))
		else:
			print("raw_R_bar incorrect. Exiting...")
			exit()
	return


def run_alpha_unit_test(input_file, community_name, list_of_qs, similarity_matrix_for_alpha_diversity, unique_species_user, unit_test=False):
	print("-------------------------------------------------------")
	print("\n\nTesting functional alpha diversity...")

	all_alpha_diversity_results = generate_final_alpha_diversity_output(input_file, community_name, input_file, list_of_qs, user_similarity_matrix=similarity_matrix_for_alpha_diversity, unique_species_user=unique_species_user)

	for i, (j1, j2) in list(all_alpha_diversity_results.items()):
		print(("\n\n\trepertoire\t%s" % i))
		for q in list_of_qs:
			if "%.3f" % j1[q] == '1.500': print(("\tq=%.1f\tpass" % q))
			else:
				print("functional_alpha incorrect. Exiting...")
				exit()
	print("\n\nConsult the unit test of recon for raw_alpha diversity\n")
	return


"""MAIN"""

if __name__ == '__main__':

	parser = ArgumentParser(description="""Calculate diversity for any qD""", formatter_class=RawTextHelpFormatter)
	pa = parser.add_argument

	# pa("-ar", "--alpha_diversity_repertoires", type=str, default=None, help="Alpha diversity will be calculated for all repertoires by default.\n\n")

	pa("-clone_distribution_in_file", "--clone_distribution_in_file", action='store_true', help='Input file is of the form clone_size tab number_of_clones of that size.\n\n') # This is the same as -c option in recon

	pa("-cn", "--community_names", type=str, default="", help="Names of communities being evaluated (>=1).\n\n")
	
	pa("-if", "--input_files", type=str, default="", help="Filenames for the community_names in the seq \t count format.\n\n")
	
	pa("-ma", "--master_filename_alpha", type=str, default="alpha_diversity_master_file.txt", help="This is the master alpha output file.\n\n")
	
	pa("-mb", "--master_filename_beta", type=str, default="beta_diversity_master_file.txt", help="This is the master beta output file.\n\n")
	
	pa("-mo", "--mode", type=str, default="beta", help="Acceptable values: alpha, beta, alpha_and_beta.\n\n")
	
	pa("-mod", "--master_output_dir", type=str, default=None, help="This is where the master output file(s) beta_diversity_master_file.txt/alpha_diversity_master_file.txt will be written.\n\n")
	
	pa("-qs", "--list_of_qs", type=str, default= "[0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., inf]", help="List of qs for which diversity is to be calculated.\n\n")
	
	pa("-rf", "--recon_files", type=str, default="", help="Filenames for the community_names in the seq \t count format. This is usually the file with all sequences.\n\n")

	pa("-u", "--unit_test", action="store_true", help="Run unit tests.\n\n")

	pa("-uq", "--unique_species_user", type=str, default=None, help="Text file that has species in the same order as used for calculating --user_similarity_matrix_file/-Z. This is the column/row header for the --user_similarity_matrix_file/-Z. Ensure that the two headers are in the same order.\nSee manual for details and example.\n\n")

	pa("-v", "--verbose", action="store_true", help="Be verbose.\n\n")
	
	pa("-weights", "--weights", type=str, default="[0.5, 0.5]", help="List of weights for r1 and r2.\n\n")
	
	pa("-Z", "--user_similarity_matrix_file", type=str, default=None, help="File containing user similarity matrix. Two options:\ni\tNumpy (.npy) matrix file\nii\tComma-separated (.csv) file with rows/columns corresponding to input matrix.\nMust also supply a file with the column/row header (note these headers must be identical; i.e., species must be in the same order in both row and column dimensions); this is supplied as -uq/--unique_species_user.\nSee manual for details and example.\n\n")

	pa("-ZF", "--user_similarity_function", type=str, default="get_similarity_matrix,%s/get_fast_similarity.py" % path_to_morty, help="Comma-separated string of length 2 in the order:\ni\tName of the python function used to calculate all-against-all similarity. This function \n\t- %s %s accept a list of species for all-against-all similarity is to be calculated\n\t- %s %s return the complete similarity matrix\nii\tPath to the .py file from which this function will be imported\nOnly one function per run can be used. See manual for details and an example.\n\n" % ( "\u0332".join("must"),  "\u0332".join("only"), "\u0332".join("must"),  "\u0332".join("only") ) )

	args = parser.parse_args()
	globals().update(vars(args))

	#			  			  #
	# *---START UNIT TEST---* #
	#			  			  #
	if unit_test:

		"""Generate two files unit_test_1.txt and unit_test_2.txt on the fly"""
		seqs_  = ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF"]
		
		with open("unit_test_1.tmp", "w") as f:
			for i in seqs_[:3]: f.write("%s\t1\n" % i)

		with open("unit_test_2.tmp", "w") as f:
			for i in seqs_[3:]: f.write("%s\t1\n" % i)

		with open("unique_seqs_user_test_file.tmp", "w") as f:
			for i in seqs_: f.write("%s\n" % i)

		with open("unique_seqs_user_test_file_1.tmp", "w") as f:
			for i in seqs_[:3]: f.write("%s\n" % i)

		with open("unique_seqs_user_test_file_2.tmp", "w") as f:
			for i in seqs_[3:]: f.write("%s\n" % i)

		input_files_list = ['unit_test_1.tmp', 'unit_test_2.tmp']
		community_names_list = ['unit_test_file_1', 'unit_test_file_2']
		list_of_qs = [0., 1., 2., 3, 3.5]
		weights = [0.5, 0.5]
		beta_diversity=True
		
		similarity_matrix_for_beta_diversity = array([ 
			[1.,  0.5, 0.5, 0.7, 0.7, 0.7], 
			[0.5, 1.,  0.5, 0.7, 0.7, 0.7], 
			[0.5, 0.5, 1.,  0.7, 0.7, 0.7], 
			[0.7, 0.7, 0.7, 1,   0.5, 0.5],
			[0.7, 0.7, 0.7, 0.5, 1,   0.5], 
			[0.7, 0.7, 0.7, 0.5, 0.5, 1  ] ])

		run_beta_unit_test(input_files_list, community_names_list, list_of_qs, similarity_matrix_for_beta_diversity, "unique_seqs_user_test_file.tmp", unit_test=True)

		similarity_matrix_for_alpha_diversity = array([
			[1.,  0.5, 0.5],
			[0.5, 1.,  0.5],
			[0.5, 0.5, 1. ] ])

		run_alpha_unit_test(input_files_list[:1], community_names_list[:1], list_of_qs, similarity_matrix_for_alpha_diversity, "unique_seqs_user_test_file_1.tmp", unit_test=True)
		run_alpha_unit_test(input_files_list[1:], community_names_list[1:], list_of_qs, similarity_matrix_for_alpha_diversity, "unique_seqs_user_test_file_2.tmp", unit_test=True)
		
		for f in ["unit_test_1.tmp", "unit_test_2.tmp", "unique_seqs_user_test_file.tmp", "unique_seqs_user_test_file_1.tmp", "unique_seqs_user_test_file_2.tmp"]:
			os.remove(f)
		
		exit()
		#			  			#
		# *---END UNIT TEST---* #
		#			  			#

	"""Start main code"""

	date_=strftime("%Y-%m-%d %H:%M:%S")
	print(""); print((strftime("%Y-%m-%d %H:%M:%S")))
	
	if input_files =="":
		print("\nProvide input\nExiting...")
		exit()
	
	start_time=time()
	run_id = ''.join(choice(ascii_letters + digits) for i in range(5))
	print(("run_id:\t%s\n" % run_id))

	if "dir_name" in argv:
		print("DEPRECATION WARNING: Command-line option --dir_name has been deprecated and will have no effect on this run")
	if "-rn" in argv or "--repertoire_names" in argv:
		print("Command line option to specify community names is now -cn/--community_names.\nExiting...")
		exit()

	"""Initialize"""
	alpha_diversity=False
	beta_diversity=False
	user_similarity_matrix=None # for when similarity matrix comes from user

	if mode=="alpha": alpha_diversity=True
	elif mode=="beta": beta_diversity=True
	elif mode=="alpha_and_beta":
		alpha_diversity=True
		beta_diversity=True
	else:
		print("Specify mode. Exiting...")
		exit()

	if unique_species_user and not user_similarity_matrix_file:
		print("Provide the similarity matrix file (--user_similarity_matrix_file/-Z) that corresponds. See help and manual for details.\nExiting...")
		exit()

	"""If the user supplies a simiarity matrix in a npy/csv file"""
	if user_similarity_matrix_file:
		print("%s: Using user-generated similarity matrix.\n" % "\u0332".join("NOTE"))
		print("User-generated similarity matrix file:\t%s" % user_similarity_matrix_file)
		print("User-generated unique species order file:\t%s" % unique_species_user)
	
		if not unique_species_user:
			print("Provide a text file (--unique_species_user/-uq) with the species in the order used to generate %s. See manual for details.\nExiting..." % user_similarity_matrix_file)
			exit()
	
		_, file_extension = os.path.splitext(user_similarity_matrix_file)

		if file_extension==".npy":
			user_similarity_matrix = np.load(user_similarity_matrix_file)
	
		elif file_extension==".csv":
			user_similarity_matrix=[]
			with open(user_similarity_matrix_file_) as f:
				for line in f:
					line=str.strip(line, "\n").split(",")
					user_similarity_matrix.append(line)
		else:
			print("-Z/--user_similarity_matrix_file should be .npy or .csv. See help and manual for details.\nExiting...")
			exit()

		"""Check if the input matrix is square"""
		try:
			assert all (len(row_) == len(user_similarity_matrix) for row_ in user_similarity_matrix)
		except AssertionError:
			print("Similarity matrix in user_similarity_matrix_file (%s) is not a square matrix. All-against-all similarity matrix should be square.\nExiting..." % user_similarity_matrix_file)
			exit()
		
		"""Make all values float"""
		user_similarity_matrix = np.array(user_similarity_matrix, dtype=np.float32)

	if not user_similarity_matrix_file:
		"""Specifcy the custom function used for similarity calculation. By default, this is fast_similarity as used in Arora and Arnaout BiorXiv 2020"""

		user_function_, function_file_path = user_similarity_function.split(",")
		print("Similarity matrix will be generated by:\nFunction:\t%s\nImported from:\t%s" % (user_function_, function_file_path))

		if len(function_file_path.split("/")) == 1:
			function_filename=function_file_path

		else:
			# sys.path.insert(0, function_file_path.split("/"))
			function_filename=function_file_path.split("/")[-1]
			sys.path.insert(0, "/".join(function_file_path.split("/")[:-1]))
		

		function_filename_, _ = os.path.splitext(function_filename)
		user_similarity_module = importlib.import_module(function_filename_)
		user_similarity_function = getattr(user_similarity_module, user_function_)	

	"""Set output files"""
	if not master_output_dir:
		print("Provide path to master output files (--master_output_dir/-mod)\nExiting...")
		exit()

	master_filename_beta = "%s/%s" % (master_output_dir, master_filename_beta)
	master_filename_alpha = "%s/%s" % (master_output_dir, master_filename_alpha)

	recon_files_list = recon_files.split(',')
	community_names_list = community_names.split(',')
	if mode=="beta" and len(community_names_list) !=2:
		print("")

	input_files_list = input_files.split(',')

	set_cython_seed(randint(0, 2147483647)) # This is cython RAND_MAX

	code_version = argv[0]
	datestamp = subprocess.getoutput("date +%m%d%Y")

	# Process list of qs
	list_of_qs = list_of_qs.replace("inf", "float('inf')")
	list_of_qs = eval(list_of_qs)

	weights = eval(weights)

	repertoire_name_to_file_name = list(zip(community_names_list, input_files_list))

	# alpha_diversity_communities_list = alpha_diversity_repertoires.split(",")

	"""Get column names for d numbers"""
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
		
	if alpha_diversity: screen_out += "\n# ALPHA DIVERSITY\n# ---------------\n# status\t%s\n# communities\t%s\n" % ( alpha_diversity, ", ".join(community_names_list) )
	if beta_diversity: screen_out += "\n# BETA DIVERSITY\n# ---------------\n# status\t%s\n# communities\t%s\n" % ( beta_diversity, ", ".join(community_names_list) )

	for i,j in repertoire_name_to_file_name:
		screen_out += "# Input for %s:\t%s\n" % (i, j)

	if verbose: print(screen_out)
	
	if alpha_diversity:
		community_ = community_names_list[0]

		all_alpha_diversity_results = generate_final_alpha_diversity_output(input_files_list, community_names_list, recon_files_list, list_of_qs, user_similarity_matrix=user_similarity_matrix, user_function=user_similarity_function, unique_species_user=unique_species_user)

		outstr_alpha = ''

		functional_alpha_diversity_results, raw_alpha_diversity_results = all_alpha_diversity_results[community_]
			
		functional_alpha_diversity_results = { str(q)+"Ds" : "%.3e"%functional_alpha_diversity_results[q] for q in list(functional_alpha_diversity_results.keys()) }
		raw_alpha_diversity_results = { str(q)+"D" : "%.3e"%raw_alpha_diversity_results[q] for q in list(raw_alpha_diversity_results.keys()) }

		outstr_alpha += "%s\talpha\t%s\t" % (run_id, community_)
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
		print("Output for run %s written to:\t%s" % (run_id, master_filename_alpha))
		
	if beta_diversity:
		if verbose: 
			print("\nCalculating beta diversity")
		
		if float('inf') in list_of_qs:
			list_of_qs_beta = [ q for q in list_of_qs if q != float("inf") ]
		elif float('-inf') in list_of_qs:
			list_of_qs_beta = [ q for q in list_of_qs if q != float("-inf") ]
		else: 
			list_of_qs_beta = list_of_qs
		
		functional_B_bar, functional_R_bar, functional_beta_bar, functional_rho_bar = generate_final_beta_diversity_output(input_files_list, community_names_list, list_of_qs_beta, user_similarity_matrix=user_similarity_matrix, user_function=user_similarity_function, unique_species_user=unique_species_user)

		functional_B_bar = { str(q)+"Ds" : "%.3e"%functional_B_bar[q] for q in list(functional_B_bar.keys()) }
		functional_R_bar = { str(q)+"Ds" : "%.3e"%functional_R_bar[q] for q in list(functional_R_bar.keys()) }

		raw_B_bar, raw_R_bar, raw_beta_bar, raw_rho_bar = generate_final_beta_diversity_output(input_files_list, community_names_list, list_of_qs_beta,diversity_type="raw", unique_species_user=unique_species_user)

		raw_B_bar = {str(q)+"D" : "%.3e"%raw_B_bar[q] for q in list(raw_B_bar.keys())}
		raw_R_bar = {str(q)+"D" : "%.3e"%raw_R_bar[q] for q in list(raw_R_bar.keys())}

		"""First do B_bar and R_bar"""
		outstr_B = ''; outstr_R = ''

		repertoire_1, repertoire_2 = community_names_list

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
		for community_ in community_names_list:
			other_repertoire = [i for i in community_names_list if i != community_][0]
			functional_beta_bar[community_] = { str(q)+"Ds" : "%.3e"%functional_beta_bar[community_][q] for q in list(functional_beta_bar[community_].keys()) }
			raw_beta_bar[community_] = { str(q)+"D" : "%.3e"%raw_beta_bar[community_][q] for q in list(raw_beta_bar[community_].keys()) }

			outstr_beta = ''
			if float('inf') in list_of_qs:
				functional_beta_bar[community_]['infDs'] = float('nan')
				raw_beta_bar[community_]['infD'] = float('nan')
			#
			if float('-inf') in list_of_qs:
				functional_beta_bar[community_]['-infDs'] = float('nan')
				raw_beta_bar[community_]['-infD']  = float('nan')

			functional_raw_beta_bar = {}
			for k_ in list(functional_beta_bar.keys()):
				functional_raw_beta_bar[k_] = dict( list(functional_beta_bar[k_].items()) + list(raw_beta_bar[k_].items()) )

			outstr_beta += "%s\tbeta_bar\t%s\t%s\t" % (run_id, community_, other_repertoire)
			outstr_beta += "%s\t" % dict(functional_raw_beta_bar[community_])
			outstr_beta += "%s\tpython %s\n" % (date_, " ".join(argv))
			outstr_beta_measures += outstr_beta

		for community_ in community_names_list:
			
			other_repertoire = [i for i in community_names_list if i != community_][0]
			
			functional_rho_bar[community_] = { str(q)+"Ds" : "%.3e"%functional_rho_bar[community_][q] for q in list(functional_rho_bar[community_].keys()) }
			raw_rho_bar[community_] = { str(q)+"D" : "%.3e"%raw_rho_bar[community_][q] for q in list(raw_rho_bar[community_].keys()) }

			outstr_rho = ''
			if float('inf') in list_of_qs:
				functional_rho_bar[community_]['infDs'] = float('nan')
				raw_rho_bar[community_]['infD']  = float('nan')
			#
			if float('-inf') in list_of_qs:
				functional_rho_bar[community_]['-infDs'] = float('nan')
				raw_rho_bar[community_]['-infD']  = float('nan')

			functional_raw_rho_bar = {}
			for k_ in list(functional_rho_bar.keys()):
				functional_raw_rho_bar[k_] = dict( list(functional_rho_bar[k_].items()) + list(raw_rho_bar[k_].items()) )

			outstr_rho += "%s\trho_bar\t%s\t%s\t" % (run_id, community_, other_repertoire)
			outstr_rho += "%s\t" % dict(functional_raw_rho_bar[community_])
			outstr_rho += "%s\tpython %s\n" % (date_, " ".join(argv))
			outstr_beta_measures += outstr_rho

		if verbose: print(("\n\n%s" % outstr_beta_measures))

		if os.path.isfile(master_filename_beta): action_ = 'a'
		else: action_ = 'w'
		with open(master_filename_beta, action_) as f: f.write(outstr_beta_measures)
		print("Output for run %s written to:\t%s" % (run_id, master_filename_beta))

	print(""); print((strftime("%Y-%m-%d %H:%M:%S")))
	if verbose: print(("# time_taken\t%.2fs" % (time() - start_time)))