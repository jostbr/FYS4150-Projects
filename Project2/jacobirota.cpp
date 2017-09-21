//Jacobi Rotation Algorithm


#include <iostream>
#include <cmath>
#include <armadillo>
#include <ctime>

using namespace std;
//using namespace arma;

void fill_array(arma::mat& A, int n){
    double rho_0 = 0.0;
    double rho_n = 5.0;
    arma::vec rho(n+1);
    rho(0) = rho_0;
    rho(n) = rho_n;

    double h_step = (rho_n - rho_0)/n;
    double hh =h_step*h_step;

    for (int i=1; i<n; i++){
        rho(i) = rho_0 + i*h_step;
    }

    arma::vec diag_el(n+1);
    for (int i=0; i<n; i++){
        diag_el(i)= (2.0/hh) + (rho(i)*rho(i));
    }

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
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
            } else {
                t = -1.0/( -tau + sqrt(1.0 + tau*tau));
            }
            c = 1/sqrt(1+t*t);
            s = c*t;
            } else {
            c = 1.0;
            s = 0.0;
        }
        double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
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
                r_ik = V(i,k);
                r_il = V(i,l);
                V(i,k) = c*r_ik - s*r_il;
                V(i,l) = c*r_il + s*r_ik;
                }

        max = 0.0;
        max_element(A, n, k, l, max);
        iterations++;
    }
    for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (A(i,j) < 1.0e-10){
                    A(i,j)= 0.0;
                }
            }
        }
    //cout << "Number of iterations was: " << iterations << endl;
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





int main(int argc, char* argv[]){
    test_max();

    string filename = argv[1];
    int n = atoi(argv[2]);


    double epsilon = 1.0e-6;

    arma::mat A = arma::zeros(n,n);
    arma::mat V = arma::eye(n,n);       //V is matrix to contain eigenvectors

    int k, l;       //Indexes for max element

    //A.print("A= ");
    //V.print("V= ");

    fill_array(A, n);

    /*
     * //Test Armadillos Eigen solver
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    //cout << "Armadillo found eigenvalues: " << eigval << endl;
    //cout << "Armadillo found eigenvectors: " << eigvec << endl;
    */



    double max = 0.0;
    max_element(A, n, k, l, max);

    clock_t start_time = clock();

    JacobiRotation(A, V, max, epsilon, n, k, l);

    clock_t end_time = clock();
    double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    cout << "Time used: " << time_used << endl;

    //A.print("A = ");            //Should contain eigenvalues along the diagonal
    //V.print("V = ");            //Should contain eigenvectors as columns
    arma::vec eig = arma::sort(A.diag());
    eig.print();
    return 0;
}

