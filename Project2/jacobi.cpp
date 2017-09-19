
# include <iostream>
# include <armadillo>

/* Function that assumes a symmetric matrix A of dimensions NxN. Then finds the
 * maximum value in the matrix excluding values along the diagonal. */
double max_non_diag(arma::mat A, int N, int* k, int* l){
    double max_val = 0.0; // Initialize to zero (checking for abs() below)

    /* Only need to loop through upper triangular part of matrix since
     * matrix is assumed symmetric. This also ensures that k < l. */
    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){      // Ensures i < j --> k < l
            if (fabs(A(i, j)) > max_val){
                max_val = A(i, j);          // Collect max val
                *k = i;                     // Store row index for max val
                *l = j;                     // Store column index for max val
            }
        }
    }

    return max_val;
}
