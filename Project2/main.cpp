
# include <iostream>
# include <armadillo>
# include "unit_tests.h"
# include "jacobi.h"

int main(){
    std::cout << "\nEXECUTING UNIT TESTS..." << std::endl
              << "============================================" << std::endl;
    TEST_get_max_non_diag();
    TEST_get_trig_values();
    TEST_jacobi_eigen();
    std::cout << "============================================" << std::endl;

    int N = 3;
    arma::mat A = arma::zeros(N, N);    // Initial matrix
    arma::mat V = arma::eye(N, N);      // To hold eigenvectors

    /* Initialize some matrix. */
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;

    A.print("\nA = ");
    jacobi_eigen(&A, &V, N);

    A.print("\nD = ");
    V.print("\nV = ");

    return 0;
}

