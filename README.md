# Floyd in Parallel

## Notes:
 1.  The input matrix is overwritten by the matrix of lengths of shortest
     paths.
 2.  Edge lengths should be nonnegative.
 3.  If there is no edge between two vertices, the length is the constant
     INFINITY.  So input edge length should be substantially less than
     this constant.
 4.  The cost of travelling from a vertex to itself is 0.  So the adjacency
     matrix has zeroes on the main diagonal.
 5.  No error checking is done on the input.
 6.  The adjacency matrix is stored as a 1-dimensional array and subscripts
     are computed using the formula:  the entry in the ith row and jth
     column is mat[i*n + j]
 7.  Use the compile flag -DSHOW_INT_MATS to print the matrix after its
     been updated with each intermediate city.
 8.  Block partition used to divide up the matrix.