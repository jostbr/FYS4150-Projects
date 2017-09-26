
# include <iostream>
# include <cmath>
# include <armadillo>

/* Function that assumes a symmetric matrix A of dimensions NxN. Then finds the
 * maximum value in the matrix excluding values along the diagonal. */
double get_max_non_diag(arma::mat A, int N, int* k, int* l){
    double max_val = 0.0; // Initialize to zero (checking for abs() below)

    /* Only need to loop through upper triangular part of matrix since
     * matrix is assumed symmetric. This also ensures that k < l. */
    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){      // Ensures i < j --> k < l
            if (fabs(A(i,j)) > fabs(max_val)){
                max_val = fabs(A(i,j));
                *k = i;                     // Store row index for max val
                *l = j;                     // Store column index for max val
            }
        }
    }

    return max_val;
}

/* Function that computes the values for cosine(theta) and sine(theta) needed for
 * the transformation matrix. The values are computed by requiring the transformed
 * element B(k,l) to be zero. This puts a constraint on theta. */
void get_trig_values(arma::mat A, int k, int l, double* cosine, double* sine){
    if (A(k,l) != 0.0){         // More work needed if A(k,l) != 0
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));              // Collecting variables
        double tangent_1 = -tau + sqrt(1 + pow(tau, 2.0));      // Largest root
        double tangent_2 = -tau - sqrt(1 + pow(tau, 2.0));      // Smallest root

        *cosine = 1/(sqrt(1 + pow(tangent_2, 2.0)));    // Trig identity for cosine
        *sine = tangent_2*(*cosine);                    // Trig identity for sine
    }

    else {      // Otherwise the values simplify
        *cosine = 1.0;
        *sine = 0.0;
    }
}

/*void transform(&A, s, c, k, l){
    // Function that performs orthogonal transformation
}*/

void jacobi_master(arma::mat A, int N){
    A.print("A = ");

    int k, l;
    double sine, cosine;

    double eps = pow(10.0, -8);
    int max_iter = (double)(N*N*N);
    double curr_max = get_max_non_diag(A, N, &k ,&l);
    get_trig_values(A, k, l, &cosine, &sine);
}
