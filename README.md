"""
Attributions:

- add_column_placeholders, parse_file_contents, ends_are_ok, extract_cdr3 and convert_fasta (with modification by Harry M. Burke) from Arnaout lab [Ramy Arnaout & Rohit Arora]
- create_recon_input_from_file, check_file_overlap, recon_pipeline, run_recon_analyses, test_resample, and separate_data_files written by Harry M. Burke as part of Arnaout lab.
- Run this script from within the same folder as recon_v2.5.py.
- You need provided style.css, plot_clone_size_distribution_ref.js (rename the file if necessary), d3.js,  error_bar_parameters.txt, and simlib.cpython-38-darwin.so on your desktop.
- Raw data folders should have no extraneous files in them aside from _recon.txt files.
- To compile simlib.pyx (cython file), run setup.py. The output compiled file you need to run is simile.cpython-38-darwin.so which is in the morty directory.

To run unit test:
morty.py -u
-> run with two files unit_test_1.txt and unit_test_2.txt within the same folder as morty.py

To do in v10:
-> There is a significant speed up to be made when calculating alpha dn beta diversity,  by doing the following:
Start calculating similarity matrix for beta diversity between two repertoires (R1 and R2; disjoint). While calculating similarity matrix, once you have reached the end of R1 check if the alpha diversity of R1 is to be calculated. If so, get the indices  of similarity matrixrelevant to R1 (i.e. ignore the similarity between R1 and R2 seqs) and throw it at the alpha diversity calculation , get the result and print it. Do the ssame for R2. If there are common seqs, book-keeping becomes slightl;y more complicated, but manageable.
-> Code for the edge case where one repertoire is a subset of another (not sure if this is needed anymore, but check)
-> Check the impact of weight ratio on beta diversity (--weights)
-> Tell user to use normalized P (i.e. P_bar) and not counts

v9:
- converted the files from python v2 to v3 

v8:
- Introduced if __name__ == '__main__' to be able to use morty as a module
- NOTE: alpha div unit test should be fixed

v7:
- Connecting alpha and beta diversities to speed up the code.

v6:
- Re-wrote almost the entire beta diversity part of the code to avoid fatal memory issues involving big arrays

v5:
- Double checked and improved the unit test

v4:
- Removed the "sample_size" option and the related code. We have found it safer to pre-sample and save it because it helps with repeatability which is what we need right now
- Added beta diversity.

v3:
- Added a --sampled_infile/-sin option which will now accept the sampled input file for diversity. We realized that naive diversity calculations and corrections on them (recon) should use the entire repertoire to get the best estimates
- We do not calculate naive diversity explicitly anymore. We get the naive diversity values from recon

v2:
- Harry added a function to read in the input file

- Implemented "no_recon" switch

  """