
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

    if (get_max_non_diag(A, N, &k, &l) == fabs(max_val)){
        std::cout << "Max non-diagonal element found ====> TEST PASSED!" << std::endl;
    }

    else {
        std::cout << "Max non-diagonal element NOT found ==== > TEST FAILED!" << std::endl;
    }
}

/* Function that acts as a unit test for get_angles() function. It uses a 3x3 matrix
 * with known (hand computed) angle for first iteration in the jacobi algorithm. Then
 * checks if get_trig_values() computes the same values for sine and cosine. */
void TEST_get_trig_values(){
    int N = 3;
    double eps = 1.0E-8;            // Some tolerance in case values aren't exact
    double k = 1, l = 2;            // Known indices for max non-diag element
    double sine, cosine;            // To hold return values from get_trig_values()

    arma::mat A = arma::zeros(N, N);
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;

    get_trig_values(A, k, l, &cosine, &sine);    // Compute trig values for transf. matrix

    if ((cosine - 1.0/sqrt(5) < eps) && (sine - (-2.0/sqrt(5)) < eps)){
        std::cout << "Correct transf. matrix found ======> TEST PASSED!" << std::endl;
    }

    else {
        std::cout << "Correct transf. matrix NOT found ======> TEST FAILED!" << std::endl;
    }
}


/* Function that acts as a unit test for function jacobi_eigen() by using a 3x3 matrix
 * with known eigenvalues and lets jacobi_eigen() compute the eigenvalues and then
 * compares with the corrects ones. Since jacobi_eigen() calls on get_max_non_diag()
 * and get_trig_values(), this test function also implicitly tests these two. */
void TEST_jacobi_eigen(){
    int N = 3;
    double eig_0 = 5.51711;     // Correct values for the eigenvalues
    double eig_1 = -0.113538;   // Correct values for the eigenvalues
    double eig_2 = 1.59642;     // Correct values for the eigenvalues
    double eps = 1.0E-5;        // Some tolerance in case values aren't exact

    arma::mat A = arma::zeros(N, N);    // To hold diagonal matrix with eigvals
    arma::mat V = arma::eye(N, N);    // Not used in this function
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;

    jacobi_eigen(&A, &V, N);    // Solve for eigenvalues

    if (A(0,0) - eig_0 > eps){
        std::cout << "Eigenvalue A(0,0) is wrong ====> TEST FAILED!" << std::endl;
    }

    else if (A(1,1) - eig_1 > eps){
        std::cout << "Eigenvalue A(1,1) is wrong ====> TEST FAILED!" << std::endl;
    }

    else if (A(2,2) - eig_2 > eps){
        std::cout << "Eigenvalue A(2,2) is wrong ====> TEST FAILED!" << std::endl;
    }

    else {
        std::cout << "Correct eigenvalues found =========> TEST PASSED!" << std::endl;
    }
}


/* Functioj that acts as a unit test in order to check if the jacobi implementation
 * preserves orthogonality through the transformations performed in jacobi_eigen(). */
void TEST_orthogonality(){
    int N = 5;
    arma::mat A = arma::zeros(N, N);    // Matrix to solve for
    arma::mat V = arma::eye(N, N);      // Matrix to hold eigenvectors

    /* Fill A with some values. */
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i,j) = i*j*j - i;    // Some arbitrary values
            }

            else {
                A(i,j) = i*j - j;   // Some arbitrary values
                A(j,i) = A(i, j);  // Make sure matrix is symmetric
            }
        }
    }

    jacobi_eigen(&A, &V, N);    // Solve eigenvalue problem

    bool test_failed = false;
    double eps = 1.0E-10;       // Tolerance for "close-to-zero"
    double inner_prod = 0.0;


    for (int c_1 = 0; c_1 < N; c_1++){      // Loop through cols in V.
        for (int c_2 = 0; c_2 < N; c_2++){  // Loop through cols in V per col in V
            inner_prod = 0.0;

            for (int k = 0; k < N; k++){
                inner_prod += V(k,c_1)*V(k,c_2);    // Compute inne product
            }

            if (inner_prod > eps && c_1 != c_2){
                test_failed = true;
            }

            if (inner_prod - 1 > eps && c_1 == c_2){
                test_failed = true;
            }
        }
    }

    if (test_failed){
        std::cout << "Orthogonality is NOT preserved ========> TEST FAILED" << std::endl;
    }

    else {
        std::cout << "Orthogonality is preserved ========> TEST PASSED!" << std::endl;
    }
}
