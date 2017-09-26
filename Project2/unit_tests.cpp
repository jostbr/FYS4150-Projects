
# include <iostream>
# include <armadillo>
# include "jacobi.h"

/* Function that acts as unit test for max_non_diag function. It uses a 5x5 matrix
 * with known largest non-diag element, then tests if max_non_diag() returns the
 * known answer. The user is the informed whether or not the test is successful. */
void TEST_get_max_non_diag(){
    int k, l, N = 5;
    double max_val = -17.4;              // Choose a max value
    arma::mat A = arma::zeros(N, N);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i, j) = i*j*j - i;    // Some arbitrary values
            }

            else {
                if (i*j == 12){
                    A(i, j) = max_val;  // Fill max value here
                }

                else {
                    A(i, j) = i*j - j;   // Some arbitrary values
                }
            }

            A(j, i) = A(i, j);  // Make sure matrix is symmetric
        }
    }

    /* Fill matrix symmetriclly with some values. */
    /*A(0,0) = 56.0; A(1,1) = 43.2; A(2,2) = 1.2;
    A(0,1) = 4.0; A(1,0) = 4.0;
    A(2,0) = max_val; A(0,2) = max_val;
    A(2,1) = 1.8; A(1,2) = 1.8;*/

    //A.print("A = ");

    if (get_max_non_diag(A, N, &k, &l) == fabs(max_val)){     // If function behaves nicely
        std::cout << "Max non-diag found ====> TEST PASSED!" << std::endl;
    }

    else {      // Otherwise notify test failed
        std::cout << "Max non-diag NOT found ==== > TEST FAILED!" << std::endl;
    }
}
