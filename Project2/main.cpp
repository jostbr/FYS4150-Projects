
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

    /*int N = 3;
    arma::mat A = arma::zeros(N, N);    // Initial matrix
    arma::mat V = arma::eye(N, N);      // To hold eigenvectors

    Initialize some matrix.
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;

    A.print("\nA = ");
    jacobi_eigen(&A, &V, N);*/

    int N = 10;
    int N_domain = N+2;
    double rho_0 = 0.0;
    double rho_max = 4.0;
    double h = (rho_max - rho_0)/(N_domain - 1);
    double rho_i, V_i;

    arma::mat A = arma::zeros(N, N);
    arma::mat V = arma::eye(N, N);

    for (int i = 0; i < N; i++){
        rho_i = rho_0 + (i+1)*h; // Only interior points
        V_i = rho_i*rho_i;

        A(i,i) = 2/(h*h) + V_i;

        if (i < N-1){
            A(i,i+1) = -1/(h*h);
        }

        if (i > 0){
            A(i,i-1) = -1/(h*h);
        }
    }

    //A.print("\nA = ");

    std::cout << rho_i << std::endl;

    jacobi_eigen(&A, &V, N);

    arma::vec eig_vec = arma::sort(A.diag());
    eig_vec.print("eig = ");
    //A.print("\nD = ");
    //V.print("\nV = ");

    return 0;
}

