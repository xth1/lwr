# LICENSE

   This software was developed by Zanoni Dias, Ulisses Dias

   It should not be redistributed or used for any commercial purpose
   without written permission from authors

   release date: nov 17, 2011

 This software is experimental in nature and is
 supplied "AS IS", without obligation by the authors to provide
 accompanying services or support.  The entire risk as to the quality
 and performance of the Software is with you. The authors
 EXPRESSLY DISCLAIM ANY AND ALL WARRANTIES REGARDING THE SOFTWARE,
 WHETHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES
 PERTAINING TO MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.


 If you use this software anytime in your work, please cite the
 following papers:

 Arruda, T. S. ; Dias, Ulisses ; Dias, Z . Heuristics for the Sorting
 by Length-Weighted Inversions Problem on Signed Permutations. In:
 First International Conference on Algorithms for Computational
 Biology, 2014, Tarragona, Spain. Algorithms for Computational
 Biology, 2014. v. 8542. p. 59-70.

 Arruda, T. S. ; Dias, Ulisses ; Dias, Z . Heuristics for the Sorting
 by Length-Weighted Inversion Problem. In: ACM Conference on
 Bioinformatics, Computational Biology and Biomedical Informatics,
 2013, Maryland, U.S.A.. Proceedings of the 4th ACM Conference on
 Bioinformatics, Computational Biology and Biomedical Informatics,
 2013. p. 499-508.


# AUTHORS

- Thiago da Silva Arruda <thiago.arruda@students.ic.unicamp.br>
- Ulisses Martins Dias   <udias@ic.unicamp.br>
- Zanoni Dias            <zanoni@ic.unicamp.br>

Last updated: November 17, 2014


# REQUIREMENTS


The program needs an installed c++ distribution. We advice the use of g++ 

We have tested the programs on linux (ubuntu), and mac OSX.


# INSTALLATION 


You should simply execute the following command inside the folder
lwr-pack.

$ make 


# USAGE 


lwr_test_metaheuristic permutation number_of_iteration inversion_limit frame_limit is_signed cost_function {-r=[[...]]}.

Parameters description:
 1. permutation: the input permutation.

 2. number_of_iteration: number of iterations that will
be performed by the Metaheuristic.

 3. inversion_limit: limit (upper bound) for the number of inversions
that will be selected at each iteration of the heuristic for building solutions.

 4. frame_limit: percentege (in the interval (0,100]) of frames that
will be selected at each iteration of the Metaheuristic.

 5. is_signed: 1 if the input permutation is signed, 0 otherwise.

 6. cost_function: 1 for linear function, 2 for quadratic function and 3 logaritmic function.

 7. -r=[[...]]: Optional parameter to give a seed solution as a input parameter


# OUTPUT  


Cost Num_Inversions [Scenario]

#EXAMPLE  


$ ./lwr 10,9,8,7,6,5,4,3,2,1 50 20 75 1 1 -r=[[1,10],[10,10],[9,9],[8,8],[7,7],[6,6],[5,5],[4,4],[3,3],[2,2],[1,1]]

$ 20 11 [[1,10],[10,10],[9,9],[8,8],[7,7],[6,6],[5,5],[4,4],[3,3],[2,2],[1,1]]

