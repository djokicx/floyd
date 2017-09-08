/* File:      floyd_parallel.c
 * 
 * Author:    Dejan Djokic
 * Section:   CS220 - 02
 * Mentor:    Peter Pacheco
 *
 * Purpose:   Implement, parallely, Floyd's algorithm for solving the
 *            all-pairs shortest path problem.
 * 
 * Compile:   mpicc -g -Wall -o fp floyd_parallel.c
 * Run:       mpiexec -n<number of processes> fp
 *
 * Input:     n, the number of vertices
 *            mat, the adjacency matrix
 * Output:    a matrix showing the costs of the shortest paths
 *
 * Notes:
 * 1.  The input matrix is overwritten by the matrix of lengths of shortest
 *     paths.
 * 2.  Edge lengths should be nonnegative.
 * 3.  If there is no edge between two vertices, the length is the constant
 *     INFINITY.  So input edge length should be substantially less than
 *     this constant.
 * 4.  The cost of travelling from a vertex to itself is 0.  So the adjacency
 *     matrix has zeroes on the main diagonal.
 * 5.  No error checking is done on the input.
 * 6.  The adjacency matrix is stored as a 1-dimensional array and subscripts
 *     are computed using the formula:  the entry in the ith row and jth
 *     column is mat[i*n + j]
 * 7.  Use the compile flag -DSHOW_INT_MATS to print the matrix after its
 *     been updated with each intermediate city.
 * 8.  Block partition used to divide up the matrix.
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

const int INFINITY = 1000000;

void Floyd(int n, int local_mat[], int my_rank, int p, MPI_Comm comm);

void Read_matrix(int n, int local_mat[], int my_rank, int p, MPI_Comm comm);
void Print_matrix(int n, int local_mat[], int my_rank, int p, MPI_Comm comm);
int min(int x, int y);

int main(int argc, char* argv[]) {
   int  n, p, my_rank;
   MPI_Comm comm;
   int* local_mat;

   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm, &my_rank);

   printf("Proc %d > This is a processor test, number %d\n", my_rank, my_rank);

   if (my_rank == 0) {
      printf("How many vertices?\n");
      scanf("%d", &n);
      printf("Enter the matrix\n");
      printf("The solution is:\n");
   }
   MPI_Bcast(&n, 1, MPI_INT, 0, comm);

   local_mat = malloc(n*n/p*sizeof(int));

   Read_matrix(n, local_mat, my_rank, p, comm);

   Floyd(n, local_mat, my_rank, p, comm);

   Print_matrix(n, local_mat, my_rank, p, comm);

   free(local_mat);
   MPI_Finalize();

   return 0;
}  /* main */

/*---------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in the adjacency matrix on process #1 and distribute
 *			  also uses MPI_Scatter
 * In args:   int n
 *            int local_mat[]
 *            int my_rank
 *            int p
 *            MPI_Comm comm
 * Out arg:   local_mat
 */
void Read_matrix(int n, int local_mat[], int my_rank, int p, MPI_Comm comm) { 
   int i, j;
   int* temp_mat = NULL;

   if (my_rank == 0) {
      temp_mat = malloc(n*n*sizeof(int));
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            scanf("%d", &temp_mat[i*n+j]);
         }  
      }
      MPI_Scatter(temp_mat, n*n/p, MPI_INT, local_mat, n*n/p, MPI_INT, 0, comm);
      free(temp_mat); //free memory space
   } 
   else {
      MPI_Scatter(temp_mat, n*n/p, MPI_INT, local_mat, n*n/p, MPI_INT, 0, comm);
   }

}  /* Read_matrix */

/*---------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix by gathering them on process #1
 * 		      also uses MPI_Gather
 * In args:   int n
 *            int local_mat[]
 *            int my_rank
 *            int p
 *            MPI_Comm comm
 */
void Print_matrix(int n, int local_mat[], int my_rank, int p, MPI_Comm comm) {
   int i, j;
   int* temp_mat = NULL;

   if (my_rank == 0) {
      temp_mat = malloc(n*n*sizeof(int));
      MPI_Gather(local_mat, n*n/p, MPI_INT, temp_mat, n*n/p, MPI_INT, 0, comm);
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            if (temp_mat[i*n+j] == INFINITY) {
               printf("i ");
            }
            else {
               printf("%d ", temp_mat[i*n+j]);
            }
         }
         printf("\n");
      }
     free(temp_mat); //free memory space
   } 
   else {
      MPI_Gather(local_mat, n*n/p, MPI_INT, temp_mat, n*n/p, MPI_INT, 0, comm);
   }
}  /* Print_matrix */

/*---------------------------------------------------------------------
 * Function:    Floyd
 * Purpose:     Apply Floyd's algorithm to the matrix mat
 * In args:     int n
 *              int local_mat[]
 *              int my_rank
 *              int p
 *              MPI_Comm comm
 * In/out arg:  local_mat:  on input, the adjacency matrix, on output
 *              lengths of the shortest paths between each pair of
 *              vertices.
 */
void Floyd(int n, int local_mat[], int my_rank, int p, MPI_Comm comm) {
   int local_int_city, int_city, j, local_city1, city2, root;
   int* row_int_city;

   row_int_city = malloc(n*sizeof(int));

   for (int_city = 0; int_city < n; int_city++) {
      root = int_city / (n/p);
      if (my_rank == root) {
         local_int_city = int_city % (n / p);
         for (j = 0; j < n; j++) {
            row_int_city[j] = local_mat[local_int_city * n + j];
         }
      }
      MPI_Bcast(row_int_city, n, MPI_INT, root, comm);
      for (local_city1 = 0; local_city1 < n/p; local_city1++) {
         for (city2 = 0; city2 < n; city2++) 
               local_mat[local_city1*n + city2] =
                  min(local_mat[local_city1*n + city2],
                     local_mat[local_city1*n + int_city]
                        + row_int_city[city2]);
      }
   # ifdef SHOW_INT_MATS
      printf("After int_city = %d\n", int_city);
      Print_matrix(mat, n);
   # endif
   }
    free(row_int_city);
}  /* Floyd */
/*---------------------------------------------------------------------
 * Function:    min
 * Purpose:     Minimum of two integers
 * In args:     int x
 *				int y
 * Out args:    the min integer
 */
 int min(int x, int y) {
 	if(x >= y) {
 		return y;
 	}
 	else{
 		return x;
 	}
 }