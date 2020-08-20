from scipy.spatial import distance
import numpy as np

def get_euclidean_distance(species_list):
	"""This function calculates all-against-all euclidean distance 
	between 2-D coordinates of two species""" 

	species_list = [ i.split("_") for i in species_list ]

	similarity_list=[]
	for species_1, x1, y1 in species_list:
		similarity_for_this_species=[]
		for species_2, x2, y2 in species_list:
			if species_1 == species_2: sim_ = 1.
			else: sim_ = distance.euclidean((float(x1),float(y1)), (float(x2), float(y2)))
			similarity_for_this_species.append(sim_)
		similarity_list.append(similarity_for_this_species)
	similarity_matrix = np.array(similarity_list)
	return  similarity_matrix

