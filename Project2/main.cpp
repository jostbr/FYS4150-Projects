
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

int main(int argc, char* argv[]){
    std::cout << "\nEXECUTING UNIT TESTS..." << std::endl
              << "============================================" << std::endl;
    TEST_get_max_non_diag();
    TEST_get_trig_values();
    TEST_jacobi_eigen();
    std::cout << "============================================" << std::endl;

    std::string filename = argv[1];
    int n = atoi(argv[2]);

    std::string fileout = filename;

    arma::mat A = arma::zeros(n,n);
    arma::mat V = arma::eye(n,n);       //V is matrix to contain eigenvectors, Orthonormal!!

    //fill_array(A, n);
    fill_array_interactive(A, n);

    //Test Armadillos Eigen solver
    clock_t start_time_ARMA = clock();

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);

    clock_t end_time_ARMA = clock();
    double time_used_ARMA = (double)(end_time_ARMA - start_time_ARMA)/CLOCKS_PER_SEC;

    std::cout << "Time used ARMADILLO: " << time_used_ARMA << std::endl;


    //cout << "Armadillo found eigenvalues: " << eigval << endl;
    //cout << "Armadillo found eigenvectors: " << eigvec << endl;

    clock_t start_time = clock();

    jacobi_eigen(&A, &V, n);

    clock_t end_time = clock();
    double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    std::cout << "Time used Jacobi: " << time_used << std::endl;

    arma::vec eig = arma::sort(A.diag());
    eig.print();

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

    std::cout << w << std::endl;
    std::cout << A(w,v) << std::endl;

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

    std::cout << w_2 << std::endl;
    std::cout << A(w_2,v_2) << std::endl;

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

    std::cout << w_3 << std::endl;
    std::cout << A(w_3,v_3) << std::endl;

    arma::vec eig_vec_1(n);
    arma::vec eig_vec_2(n);
    arma::vec eig_vec_3(n);

    //Defining the wavefunction from the eigenvectors
    for (int j=0; j<n; j++){
        eig_vec_1(j) = V(j,w)*V(j,w);
        eig_vec_2(j) = V(j,w_2)*V(j,w_2);
        eig_vec_3(j) = V(j,w_3)*V(j,w_3);
    }


    //write_results_to_file(fileout, eig, V, n);
    //write_results_to_file_plot(fileout, eig, eig_vec_1, eig_vec_2, eig_vec_3, n);
    write_results_to_file_plot(fileout, eig_vec_1, eig_vec_2, eig_vec_3, n);


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

