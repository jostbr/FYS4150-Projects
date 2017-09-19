
# include <iostream>
# include <armadillo>
# include "jacobi.h"

/* Function that acts as unit test for max_non_diag function. It uses a 3x3 matrix
 * with known largest non-diag element, then tests if max_non_diag() returns the
 * known answer. The user is the informed whether or not the test is successful. */
void test_max_non_diag(){
    int k, l, N = 3;
    double max_val = 17.0;
    arma::mat A = arma::zeros(N, N);

    /* Fill matrix symmetriclly with some values. */
    A(0,0) = 56.0; A(1,1) = 43.2; A(2,2) = 1.2;
    A(0,1) = 4.0; A(1,0) = 4.0;
    A(2,0) = max_val; A(0,2) = max_val;
    A(2,1) = 1.8; A(1,2) = 1.8;
    //A.print("A = ");

    if (max_non_diag(A, N, &k, &l) == max_val){     // If function behaves nicely
        std::cout << "Max non-diag found ====> TEST PASSED!" << std::endl;
    }

    else {      // Otherwise notify test failed
        std::cout << "Max non-diag NOT element ==== > TEST FAILED!" << std::endl;
    }
}
