
# include <iostream>
# include <armadillo>
# include <cmath>

/* Function that assumes a symmetric matrix A of dimensions NxN. Then finds the
 * maximum value in the matrix excluding values along the diagonal. */
double get_max_non_diag(arma::mat A, int N, int* k, int* l){
    double max_val = 0.0; // Initialize to zero (checking for abs() below)

    /* Only need to loop through upper triangular part of matrix since
     * matrix is assumed symmetric. This also ensures that k < l. */
    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){      // Ensures i < j --> k < l
            if (fabs(A(i, j)) > fabs(max_val)){
                max_val = fabs(A(i, j));          // Collect max val
                *k = i;                     // Store row index for max val
                *l = j;                     // Store column index for max val
            }
        }
    }

    return max_val;
}

/*void get_angles(A, k, l, &s, &c){
    // Function that computes new sin and cos for transformation matrix
}

void transform(&A, s, c, k, l){
    // Function that performs orthogonal transformation
}*/

void jacobi_master(arma::mat A, int N){
    A.print("A = ");

    int k, l;
    double eps = pow(10.0, -8);
    int max_iter = (double)(N*N*N);
    double curr_max = get_max_non_diag(A, N, &k ,&l);
}
