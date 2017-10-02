
# include <iostream>
# include <armadillo>
# include <iomanip>
# include <fstream>
# include <string>
# include <ctime>
# include "initialize.h"
# include "unit_tests.h"
# include "jacobi.h"

void write_results_to_file_plot(std::string fileout, arma::vec eig_vec_1, arma::vec eig_vec_2, arma::vec eig_vec_3, int n);

/* Main function initially runs some unit tests to verify that everything works as it should.
 * Then sets up matrices and calls on functions in initialize.cpp to initialize the matrices
 * either for the single electron case or the interacting case. Then it calls on jacobi_eigen()
 * in jacobi.cpp in order to solve for eigenvalues and eigenvectors. At last eigenvectors are
 * extracted (for corresponding eigenvalues) and then are written to file. */
int main(int argc, char* argv[]){
    if (argc != 3){
        std::cout << "BAD USAGE OF COMMAND LINE ARGUMENTS! NEED 2." << std::endl;
        std::cout << "CMD ARGS: <output-filename> <N>" << std::endl;
        std::cout << "\nTerminating program.." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "\nEXECUTING UNIT TESTS..." << std::endl
              << "==================================================" << std::endl;
    TEST_get_max_non_diag();
    TEST_get_trig_values();
    TEST_jacobi_eigen();
    TEST_orthogonality();
    std::cout << "==================================================\n" << std::endl;

    std::string filename = argv[1];
    int n = atoi(argv[2]);

    std::string fileout = filename;

    arma::mat A = arma::zeros(n,n);
    arma::mat V = arma::eye(n,n);       //V is matrix to contain eigenvectors, Orthonormal!!

    fill_array(A, n);                   // Initialize non-interacting case
    //fill_array_interactive(A, n);     // Initialize interacting case

    // Time Armadillos Eigen solver.
    arma::vec eigval;
    arma::mat eigvec;
    clock_t start_time_ARMA = clock();
    arma::eig_sym(eigval, eigvec, A);
    clock_t end_time_ARMA = clock();
    double time_used_ARMA = (double)(end_time_ARMA - start_time_ARMA)/CLOCKS_PER_SEC;

    // Time jacobi implementation.
    clock_t start_time = clock();
    jacobi_eigen(&A, &V, n);        // Solve eigenvalue problem
    clock_t end_time = clock();
    double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;

    std::cout << "Timing for " << n << "x" << n << " matrix:" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Time used ARMADILLO: " << time_used_ARMA << std::endl;
    std::cout << "Time used Jacobi: " << time_used << std::endl;

    arma::vec eig = arma::sort(A.diag());   // Sort eigenvalues in increasing order
    eig.print("\nEigenvalues= ");           // Recomend only print for small n

    //Find index of three first wavefunc
    double min_eigval = 10.0e4;
    int w, v;

    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            if(i==j){
                if (fabs(A(i,j)) < min_eigval){
                    min_eigval = fabs(A(i,j));
                    w = i;
                    v = j;
                }
            }
        }
    }

    double min_eigval_2 = 10.0e4;
    int w_2, v_2;

    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            if(i==j){
                if (fabs(A(i,j)) > min_eigval && fabs(A(i,j)) < min_eigval_2 ){
                    min_eigval_2 = fabs(A(i,j));
                    w_2 = i;
                    v_2 = j;
                }
            }
        }
    }

    double min_eigval_3 = 10.0e4;
    int w_3, v_3;

    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            if(i==j){
                if (fabs(A(i,j)) > min_eigval_2 && fabs(A(i,j)) < min_eigval_3 ){
                    min_eigval_3 = fabs(A(i,j));
                    w_3 = i;
                    v_3 = j;
                }
            }
        }
    }

    arma::vec eig_vec_1(n);
    arma::vec eig_vec_2(n);
    arma::vec eig_vec_3(n);

    //Defining the wavefunction from the eigenvectors
    for (int j=0; j<n; j++){
        eig_vec_1(j) = V(j,w)*V(j,w);
        eig_vec_2(j) = V(j,w_2)*V(j,w_2);
        eig_vec_3(j) = V(j,w_3)*V(j,w_3);
    }

    write_results_to_file_plot(fileout, eig_vec_1, eig_vec_2, eig_vec_3, n);

    std::cout << "\nThree first wavefunctions written to " << fileout << std::endl << std::endl;

    return 0;
}

void write_results_to_file_plot(std::string fileout, arma::vec eig_vec_1, arma::vec eig_vec_2, arma::vec eig_vec_3, int n){
    std::ofstream ofile;    // File object for output file
    ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "      Wavevector1:        Wavevector2:           Wavevector3:" << std::endl;
    for (int j = 0; j<n; j++){
        ofile << std::setw(20) << std::setprecision(8) << eig_vec_1(j);
        ofile << std::setw(20) << std::setprecision(8) << eig_vec_2(j);
        ofile << std::setw(20) << std::setprecision(8) << eig_vec_3(j) << std::endl;
     }

    ofile.close();
}

