

Prerequisites:

Installation of the dlib optimization routines http://dlib.net/optimization.html



Compile:

path_to_dlib="path_to_dlib_library_root_directory"

g++ -std=c++0x -O3 -Wall -I$path_to_dlib -o relcoas relcoas.cpp









##########
				Example
##########



This example is a simulated data set with L=19839 sites and 48 individuals. Individual sequences are mated according to pedigree f2.1 from Theunert et al. 2017. The data are simulated randomly and do not represent any specific population. They only serve to demonstrate the data format and how to run the software.

File 'alleles.txt' contains the derived ('1') and ancestral ('0') alleles for each of the sites and each of the individuals in an ms (Hudson 2002) like format (i.e. one individual per row, one genomic site per column). 

For example 3 individuals and 4 sites would look like this:
0110
1100
1110

The file 'metadata.txt' contains 3 columns with site_number, qk, fk separated by tabs. Note that the site number is not an acutal genomic position but only a consecutive number starting at 1. qk and fk are the site specific derived allele frequencies from the ancient samples (qk) and the reference panel population (fk).


These are the simulated relatedness coefficients for pairs of individuals. All other pairs are unrelated.

#####     Simulated relatedness r_ij for pairs of individuals i and j
i/j     1       46      47      48
1       -       0.5    0.25    0.25
46      -       -      0.5     0.5
47      -       -       -      1.0
48      -       -       -       -

The individuals were contaminated with the folliwing C rates:
0.16,0.02,0.16,0.16,0.13,0.07,0.1,0.13,0.01,0.1,0.0,0.16,0.14,0.15,0.15,0.07,0.08,0.09,0.09,0.01,0.17,0.05,0.09,0.06,0.19,0.12,0.13,0.17,0.04,0.15,0.19,0.07,0.06,0.01,0.18,,0.17,0.2,0.09,0.05,0.1,0.09,0.05,0.11,0.09,0.08,0.1,0.07


##### relcoas parameter:

-indiv:
	- a single number of N individuals if you want to calculate all possible pairs for N individuals (e.g. "-indiv 10" if N=10 in your dataset)
	- a comma separated list of numbers, indicating which individuals you want to look at (e.g. "-indiv 1,3,6,8" if you only want to look at these 4 individuals)					

-in_file:
	- the name of the input allele file in ms like format (e.g. "-in_file alleles.txt")
	
-meta_file:
	- the name of the input metadata file with qk and fk (e.g. "-meta_file metadata.txt")

-max_e:
	- upper bound value of parameter e (sequencing error) that should be  used for optimization (e.g. "-max_e 0.001")

-max_c:
	- upper bound value of parameter C (contamination rate) that should be  used for optimization (e.g. "-max_c 0.25")
	

	
##### Run relcoas with only a subset of all individuals:

./relcoas -indiv 1,2,46,47,48 -in_file alleles.txt -meta_file metadata.txt  -max_e 0.001 -max_c 0.2


##### Output:

Estimated relatedness coefficients r_ij
        1       2       46      47      48
1       -       0.00    0.42    0.16    0.17
2       -       -       0.43    0.10    0.12
46      -       -       -       0.39    0.39
47      -       -       -       -       0.95
48      -       -       -       -       -

Contamination Estimates C_i
0.01 0.15 0.11 0.15 0.13 0.05 0.04 0.10 0.04 0.06 0.14 0.12 0.16 0.05 0.19 0.06 0.14 0.20 0.07 0.17 0.14 0.05 0.18 0.05 0.10 0.14 0.10 0.11 0.20 0.05 0.13 0.13 0.20 0.15 0.19 0.14 0.14 0.02 0.19 0.16 0.09 0.02 0.08 0.14 0.13 0.07 0.09 0.05 

Sequencing Error e 
0.0001

