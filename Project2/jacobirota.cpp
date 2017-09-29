//Jacobi Rotation Algorithm


#include <iostream>
#include <cmath>
#include <armadillo>
#include <ctime>
#include <fstream>
#include <string>
#include <iomanip>


using namespace std;
//using namespace arma;

void fill_array(arma::mat& A, int n){
    double rho_0 = 0.0;
    double rho_n = 4.0;
    arma::vec rho(n+2);
    rho(0) = rho_0;
    rho(n+1) = rho_n;

    double h_step = (rho_n - rho_0)/(n+1);
    double hh =h_step*h_step;

    //cout << "h_step= " << h_step << endl;
    //cout << "hh = " << hh << endl;

    for (int i=1; i<n+1; i++){
        rho(i) = rho_0 + i*h_step;
    }

    //cout << "rho" << rho << endl;

    arma::vec diag_el(n+1);
    for (int i=0; i<n+1; i++){
        diag_el(i)= (2.0/hh) + (rho(i)*rho(i));
    }

    double off_const = -1.0/hh;

    //I am not including rho_0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j){A(i,j)=diag_el(i+1);}
            if (fabs(i-j) == 1){A(i,j)=off_const;}

        }
    }
    //diag_el.print("Diag element");
}

void fill_array_interactive(arma::mat& A, int n){
    double rho_0 = 0.0;
    //Rho max scales with the frequency omega
    double rho_n = 7.0;
    arma::vec rho(n+2);
    rho(0) = rho_0;
    rho(n+1) = rho_n;

    double h_step = (rho_n - rho_0)/(n+1);

    cout << "h_step= " << h_step << endl;

    double hh =h_step*h_step;

    cout << "hh = " << hh << endl;

    for (int i=1; i<n+1; i++){
        rho(i) = rho_0 + i*h_step;
    }
    //cout << "rho= " << rho << endl;

    //Defining the frequency, 0.01, 0.5, 1, 5
    double omega = 0.500;
    double omega_squared = omega*omega;


    arma::vec diag_el(n);

    for (int i=0; i<n; i++){
        diag_el(i)= (2.0/hh) + (rho(i+1)*rho(i+1))*omega_squared + (1.0/rho(i+1));
    }

    //diag_el.print("Diag element = ");

    double off_const = -1.0/hh;

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j){A(i,j)=diag_el(i);}
            if (fabs(i-j) == 1){A(i,j)=off_const;}

        }
    }

    //A.print("A= ");

}

void max_element(arma::mat& A, int n, int& k, int& l, double& max){
    //double max = 0.0;

    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            if (i != j){
                if (fabs(A(i,j)) > max){
                    max = (fabs(A(i,j)));
                k=i;
                l=j;
                //cout << "Max = " << max << endl;
                //cout << "k = " << k << endl;
                //cout << "l = " << l << endl;
                }
            }
            else{max = max;}
        }
    }
    //cout << "Max = " << max << endl;
    //cout << "k = " << k << endl;
    //cout << "l = " << l << endl;
}


void max_element_tridiag(arma::mat& A, int n, int& k, int& l, double& max){
    //Want to take advantage of the tridiagonal matrix form when searching for max element
    //BUT when we rotate - can off diag element become non zero?
    //YES!! this do not work!!

    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            if (fabs(i-j) == 1){
                if (fabs(A(i,j)) > max){
                    max = (fabs(A(i,j)));
                k=i;
                l=j;
                }
            }

        }
    }
    //cout << "Max = " << max << endl;
    //cout << "k = " << k << endl;
    //cout << "l = " << l << endl;
}




void JacobiRotation(arma::mat& A, arma::mat& V, double& max, double& epsilon, int n, int& k, int& l){
    int iterations = 0.0;
    int max_iterations = n*n*n;

    double s, c;
    double t, tau;

    while ( fabs(max) > epsilon && (double) iterations < max_iterations ){
        //cout << "Max value 1=" << max << endl;
        if ( A(k,l) != 0.0 ) {
            tau = (A(l,l) - A(k,k))/(2*A(k,l));
            if ( tau > 0 ) {
                //t = -tau + sqrt(1.0 + tau*tau);
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
            } else {
                //t = -tau - sqrt(1.0 + tau*tau);
                t = -1.0/( -tau + sqrt(1.0 + tau*tau));
            }
            c = 1.0/sqrt(1.0+t*t);
            s = c*t;
            } else {
            c = 1.0;
            s = 0.0;
        }
        //cout << "t=" << t <<endl;
        double a_kk, a_ll, a_ik, a_il;
        a_kk = A(k,k);
        a_ll = A(l,l);
        //A.print("A= ");
        // changing the matrix elements with indices k and l
        A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        A(k,l) = 0.0; // hard-coding of the zeros
        A(l,k) = 0.0;
        //A.print("A= ");
        // and then we change the remaining elements
        for ( int i = 0; i < n; i++ ) {
            if ( i != k && i != l ) {
                a_ik = A(i,k);
                a_il = A(i,l);
                A(i,k)= c*a_ik - s*a_il;
                A(k,i) = A(i,k);
                A(i,l) = c*a_il + s*a_ik;
                A(l,i) = A(i,l);
                }
            // Finally, we compute the new eigenvectors
                double v_ik = V(i,k);
                double v_il = V(i,l);
                //V(i,k) = c*v_ik - s*(v_il + tau*v_ik);
                //V(i,l) = c*v_il + s*(v_ik - tau*v_il);
                V(i,k) = c*v_ik - s*v_il;
                V(i,l) = c*v_il + s*v_ik;
                }

        max = 0.0;
        max_element(A, n, k, l, max);
        //max_element_tridiag(A, n, k, l, max);
        //cout << iterations << endl;
        iterations++;

    }

    //Forces all elements to be 0.0 below a treshold,
    double equal_zero = 10.0e-9;
    for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (i != j && A(i,j) < equal_zero){
                        A(i,j)= 0.0;
                    }
            }
        }

    cout << "Number of iterations was: " << iterations << endl;

}


/*UNIT TEST TO CHECK IF RIGHT MAX OFF DIAGONAL ELEMENT IS RETURNED
 * In this test the function max_element should return max element off-diagonal in a
 * known matrix B.
 * B(0,2) should be return, giving index k=0, l=2
 * If test is passed, a positiv comment is prompt in the terminal window
 */

void test_max(){
    int N = 3;
    arma::mat B(N,N);
    B(0,0) = 3.0;
    B(0,1) = 1.0;
    B(0,2) = 5.0;
    B(1,0) = 1.0;
    B(1,1) = 3.0;
    B(1,2) = -1.0;
    B(2,0) = -1.0;
    B(2,1) = -1.0;
    B(2,2) = 5.0;

    int k,l;
    double max = 0.0;

    max_element(B,N, k,l,max);
    //cout << "k= " << k << endl;
    if (k == 0 && l==2){
        cout << "Max element test passed" << endl;
    }
}

void write_results_to_file(string fileout, arma::vec eig, arma::mat V, int n){
    ofstream ofile;    // File object for output file
    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       lamda:             eigenvectors:          " << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j<n; j++){
            ofile << setw(20) << setprecision(8) << eig(i);
            ofile << setw(20) << setprecision(8) << V(i,j) << endl;
        }

    }

    ofile.close();
}

void write_results_to_file_plot(string fileout, arma::vec eig, arma::vec eig_vec_1, arma::vec eig_vec_2, arma::vec eig_vec_3, int n){
    ofstream ofile;    // File object for output file
    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "      Wavevector1:        Wavevector2:           Wavevector3:" << endl;
    for (int j = 0; j<n; j++){
        ofile << setw(20) << setprecision(8) << eig_vec_1(j);
        ofile << setw(20) << setprecision(8) << eig_vec_2(j);
        ofile << setw(20) << setprecision(8) << eig_vec_3(j) << endl;

     }



    ofile.close();
}



int main(int argc, char* argv[]){
    test_max();

    string filename = argv[1];
    int n = atoi(argv[2]);

    string fileout = filename;


    double epsilon = 10.0e-10;

    arma::mat A = arma::zeros(n,n);
    arma::mat V = arma::eye(n,n);       //V is matrix to contain eigenvectors, Orthonormal!!


    int k, l;       //Indexes for max element

    //A.print("A= ");
    //V.print("V= ");

    fill_array(A, n);
    //fill_array_interactive(A, n);


    /* Just a test array with known eig_val and eig_vec
    A(0,0) = 2.0; A(1,1) = 4.0; A(2,2) = 1.0;
    A(0,1) = 1.0; A(1,0) = 1.0;
    A(0,2) = -1.0; A(2,0) = -1.0;
    A(1,2) = -2.0; A(2,1) = -2.0;
    */



    //Test Armadillos Eigen solver
    //clock_t start_time = clock();
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);

    //clock_t end_time = clock();
    //double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;

    //cout << "Armadillo found eigenvalues: " << eigval << endl;
    //cout << "Armadillo found eigenvectors: " << eigvec << endl;


    double max = 0.0;
    max_element(A, n, k, l, max);
    //max_element_tridiag(A, n, k, l, max);

    clock_t start_time = clock();

    JacobiRotation(A, V, max, epsilon, n, k, l);

    clock_t end_time = clock();
    double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    cout << "Time used: " << time_used << endl;

    //A.print("A = ");            //Should contain eigenvalues along the diagonal
    //V.print("V = ");            //Should contain eigenvectors as columns

    arma::vec eig = arma::sort(A.diag());
    //eig.print();

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

    cout << w << endl;
    cout << A(w,v) << endl;

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

    cout << w_2 << endl;
    cout << A(w_2,v_2) << endl;

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

    cout << w_3 << endl;
    cout << A(w_3,v_3) << endl;

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
    write_results_to_file_plot(fileout, eig, eig_vec_1, eig_vec_2, eig_vec_3, n);
    delete [] & A, delete [] & V, delete[] & eig;

    return 0;
}


