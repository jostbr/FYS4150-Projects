
# include <iostream>
# include <cmath>
# include <armadillo>
# include "jacobi.h"

/* Function that acts as unit test for max_non_diag() function. It uses a 5x5 matrix
 * with known largest non-diag element, then tests if max_non_diag() returns the
 * known answer. The user is the informed whether or not the test is successful. */
void TEST_get_max_non_diag(){
    int k, l, N = 5;
    double max_val = -17.4;              // Choose a max value
    arma::mat A = arma::zeros(N, N);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i,j) = i*j*j - i;    // Some arbitrary values
            }

            else {
                if (i*j == 12){
                    A(i,j) = max_val;  // Fill max value here
                }

                else {
                    A(i,j) = i*j - j;   // Some arbitrary values
                }
            }

            A(j,i) = A(i, j);  // Make sure matrix is symmetric
        }
    }

    //A.print("A = ");

    if (get_max_non_diag(A, N, &k, &l) == fabs(max_val)){     // If function behaves nicely
        std::cout << "Max non-diag found ====> TEST PASSED!" << std::endl;
    }

    else {      // Otherwise notify test failed
        std::cout << "Max non-diag NOT found ==== > TEST FAILED!" << std::endl;
    }
}

/* Function that acts as a unit test for get_angles() function. It uses a 3x3 matrix
 * with known (hand computed) angle for first iteration in the jacobi algorithm. Then
 * checks if get_trig_values() computes the same values for sine and cosine. */
void TEST_get_trig_values(){
    int N = 3;
    arma::mat A = arma::zeros(N, N);
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;

    double eps = pow(10.0, -8.0);   // Some tolerance in case values aren't exact
    double k = 1, l = 2;            // Known indices for max non-diag element
    double sine, cosine;            // To hold return values from get_trig_values()

    get_trig_values(A, k, l, &cosine, &sine);    // Compute angles

    if ((cosine - 1.0/sqrt(5) < eps) && (sine - (-2.0/sqrt(5)) < eps)){     // If function behaves nicely
        std::cout << "Correct angles found ====> TEST PASSED!" << std::endl;
    }

    else {      // Otherwise notify test failed
        std::cout << "Correct angles NOT found ==== > TEST FAILED!" << std::endl;
    }
}
