//Jacobi Rotation Algorithm


#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

//int argc, char* argv[]
//int n = atoi(argv[1]);
int main()
{
    int n = 3;
    mat A(n,n);
    /*A(0,0) = 3.0;
    A(0,1) = 1.0;
    A(0,2) = -1.0;
    A(1,0) = 1.0;
    A(1,1) = 3.0;
    A(1,2) = -1.0;
    A(2,0) = -1.0;
    A(2,1) = -1.0;
    A(2,2) = 5.0;*/
    //R should return 2,3,6 along the diagonal

    A(0,0) = 3.0;
    A(0,1) = -1.0;
    A(0,2) = 0.0;
    A(1,0) = -1.0;
    A(1,1) = 2.0;
    A(1,2) = -1.0;
    A(2,0) = 0.0;
    A(2,1) = -1.0;
    A(2,2) = 3.0;
    A.print("A= ");
    //R should return 1,4,3 along the diagonal

    mat R = eye(n,n);
    R.print("R= ");

    int k, l;
    double max = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            if (fabs(A(i,j)) > max){
                max = fabs(A(i,j));
                k = i;
                l = j;
            }


        }
    }
    cout << "Max value =" << max << endl;
    cout << "k=" << k << endl;
    cout << "l=" << l << endl;

    double epsilon = 1.0e-8;
    double max_number_iterations =  3*n;
    int iterations = 0;
    double s, c;
    double t, tau;

    while ( fabs(max) > epsilon && (double) iterations < max_number_iterations ){
        cout << "Max value 1=" << max << endl;
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
        A.print("A= ");
        // changing the matrix elements with indices k and l
        A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
        A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
        A(k,l) = 0.0; // hard-coding of the zeros
        A(l,k) = 0.0;
        A.print("A= ");
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
                r_ik = R(i,k);
                r_il = R(i,l);
                R(i,k) = c*r_ik - s*r_il;
                R(i,l) = c*r_il + s*r_ik;
                }
        A.print("A= ");
        R.print("R= ");
        max = 10e-9;
        cout << "Max value 2=" << max << endl;
        for (int i = 0; i < n; i++){
            for (int j = i+1; j < n; j++){
                if (fabs(A(i,j)) > max){
                    max = fabs(A(i,j));
                    k = i;
                    l = j;
                    cout << "k=" << k << endl;
                    cout << "l=" << l << endl;
                    cout << "Max value 3=" << max << endl;
                    cout << "Number of iterations=" << iterations << endl;

                }
                else{
                    cout << "Jacobi rotational algorithm converged" << endl;
                }
            }
        }
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;
    R.print("R= ");
    A.print("A= ");


    return 0;

}
