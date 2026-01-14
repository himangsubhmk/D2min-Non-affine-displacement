code: D2min.c

compilation(need gsl library):  gcc D2min.c -lgsl -lgslcblas -lm

execution: ./a.out N (iframe-1) iframe 

See run_D2min.sh

This code need to read Dump-file-series [a] (See the code where it reads)
Output will be saved in a folder with filename D2min_frame”n”.dat where D2min of each particles are listed for particle id 1-N
