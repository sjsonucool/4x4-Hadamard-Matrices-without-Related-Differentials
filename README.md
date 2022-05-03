# 4x4-Hadamard-Matrices-without-Related-Differentials

All programs can be compiled with 1) cc -O3 [program_filename] or 2) cc -O3 -o [executable_name] [program_filename] 

searchRelatedDiffsandRelations.c
================================
The program generates random 4x4 Hadamard MDS matrices over GF(256) and searches for related differentials and the underlying relations satisfied by the matrix elements.

filterNoRelMatrices.c
=====================
compiling and running ./a.out i (where i can be {1,2,....,255}), the program generates file containing all 4x4 Hadamard MDS matrices M with M[0,0] = i where the elements satisfy none of the 28 relations. Let us denote the list as L_i. The size of each generated file can be approximately 250 MB and each consists of 14229600 matrices. The repository already has the list of matrices with M[0,0]=1, i.e., L_1 which can be found in the folder named Files.

verRelDiff.c
============
compiling and running ./a.out i, the program searches related differentials for each matrices generated in files by program above. In order to verify correctness of the results in the paper, it is sufficient to verify the file with M[0,0] = 1, i.e., L_1. Because all other lists L_2, L_3,...,L_255 contains matrices which are scalar multiples of Matrices in L_1.

Note:
====
One can run multiple instances of a.out in parallel for faster verification. The file with M[0,0] = 1 can be split further into subfiles and verified in parallel too!
