
# include <iostream>
# include <armadillo>
# include "jacobi.h"

double max_non_diag(arma::mat, int, int*, int*);
void test_max_non_diag();

int main(){
    //double val = 10.0;
    //jacobi_method(val);

    test_max_non_diag();

    return 0;
}

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

/* Function that acts as unit test for max_non_diag function. It uses a 3x3 matrix
 * with known largest non-diag element, then tests if max_non_diag() returns the
 * known answer. The user is the informed whether or not the test is successful. */
void test_max_non_diag(){
    int k, l, N = 3;
    double max_val = 17.0;
    arma::mat A = arma::zeros(N, N);

    A(0,0) = 56.0; A(1,1) = 43.2; A(2,2) = 1.2;
    A(0,1) = 4.0; A(1,0) = 4.0;
    A(2,0) = max_val; A(0,2) = max_val;
    A(2,1) = 1.8; A(1,2) = 1.8;
    //A.print("A = ");

    if (max_non_diag(A, N, &k, &l) == max_val-1){     // If function behaves nicely
        std::cout << "Max non-diag found ====> TEST PASSED!" << std::endl;
    }

    else {      // Otherwise notify test failed
        std::cout << "Max non-diag NOT element ==== > TEST FAILED!" << std::endl;
    }
}

