"""
To compile:

1. make sure that setup.py reads in the current version (here simlib.pyx) and has calculate as its target (i.e., without versions)

2. python setup.py build_ext --inplace
"""

import cython

from libc.stdlib cimport RAND_MAX
from libc.time cimport time
from libc.math cimport log

cdef extern from "stdlib.h":
	double drand48()
	void srand48(long int seedval)

cpdef set_cython_seed(unsigned int random_number_from_python):
	srand48( random_number_from_python )
	return

set_cython_seed(RAND_MAX)

cpdef stochastic_similarity(unsigned int x=1, mu=0.53):
	# x is the edit distance between two sequences; defaults for convenience to 1
	# mu is the mean log10(fold effect)
	# c = 1./(10**mu) = 0.1**mu
	# Sensitivity analysis should cover the range 0.45 to 0.69. See log.rtfd 091318RA "Best k"
	# If stochasticity_flag is on, then sample from an exponential probability-density function (k*exp(-k*t) ) that has mu as the average value. The average of an exponential PDF is 1/k, where k is the rate parameter; so k = 1/mu. This is why you multiply the log() by mu; it's the same as log()/k.
	#
	# We find, empirically, that stochastic sampling from either an exponential fit of the data or the data itself (using bisect()) leads to a relationship between edit distance and similarity that follows a s=kappa**x relationship (over many orders of magnitude). For example, for our original measure of s=0.1**(x*4/(mean_length)), we end up with a kappa of 0.6 in a simulation with mean CDR3 length in the repertoire = 18, with s.d. 2. I.e., s=0.6**x. In this case, it's easy to see that this will be the case for mean_length=18. This kappa differs from k (=1/mu), which is the rate constant of exponential fits to the data. The point is, even when you stochastically sample using mu, you end up with a relationship kappa**x. So we can computationally take the shortcut of simply (deterministically) using s=kappa*x.  That's the "else" below. So if stochasticity_flag == 1: use mu. We currently see no reason to do this other than to confirm that we get the same result with the appropriately chosen kappa. Otherwise (if stochasticity_flag = 0, the default), just use the passed kappa and do the faster calculation. See log.rtfd 091718RA.
	cdef:
		double random
		double s
		double c
	s = 1. # initialize similarity
	for i from 1 <= i <= x:
		random = drand48()
		c = 10**( log(random) * mu )
		s = s * c
	return s


cpdef fast_similarity(unsigned int x=1, double cost=0.51):
		return cost**x
