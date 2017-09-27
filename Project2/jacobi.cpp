#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;
ofstream ofile;

// Function to find the maximum off-diagonal element of a matrix A and the indices of this element.
void max_element(double ** A, int* k, int* l, double* max_elem, int n){
   *k = *l = 0;
   *max_elem =0.0;     // Sets initial values to zero
    for (int i=0; i<n; i++){
       for (int j = i+1; j<n; j++){

      double a = A[i][j];
          double aa = a*a;
          if (aa > *max_elem){
             *max_elem = aa;

             *k = i;
             *l = j;
           }
        }
     }
}

void func(int k, int l, int n, double ** A, double ** B, double **S){
   double t; double s; double c;
   if (A[k][l] != 0){
     double tau = (A[l][l]-A[k][k])/(2.0*A[k][l]);
     if (tau > 0){ t=1.0/(tau+sqrt(1.0+tau*tau));}
     else { t=-1.0/(-tau+sqrt(1.0+tau*tau));}
     c = 1.0/sqrt(1+t*t);
     s = t*c;
     }
   else {
      c = 1.0;
      s = 0.0;
   }
   double cc = c*c;
   double ss = s*s;
   double cs = c*s;

   for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++){
         if ( i !=k && i!=l){
            if (j==i){B[i][j] = A[i][j];}
            else if (j==k){B[i][j]=B[j][i]=A[i][k]*c-A[i][l]*s;}
            else if (j==l){B[i][j]=B[j][i]=A[i][l]*c+A[i][k]*s;}
         else {B[i][j] = A[i][j];}
         }
       }
   }
   B[k][k] = A[k][k]*cc-2*A[k][l]*cs+A[l][l]*ss;
   B[l][l] = A[l][l]*cc+2*A[k][l]*cs+A[k][k]*ss;
   B[k][l] = (A[k][k]-A[l][l])*cs+A[k][l]*(cc-ss);
   B[l][k] = -B[k][l];

   // Setting up new eigenvectors

   for (int i =0; i<n; i++){
      double s_ik = S[i][k];
      double s_il = S[i][l];
      S[i][k] = c*s_ik-s*s_il;
      S[i][l] = c*s_il+s*s_ik;

   }


}

void unit_test(double ** B, double ** S){
   // Function to check if the Jacobi method gives the right eigenvalues for a 4x4 matrix with known eigenvalues
   int n = 4;
   //Defining matrix
   double ** C = new double * [n];
   for (int i = 0; i<4; i++){
      C[i] = new double [n];}
   C[0][0] = C[3][3] = 3; C[0][1] = C[1][0] = 1; C[0][2] = C[2][0] = -1; C[0][3]=C[3][0] = 2; C[1][1] = 4; C[2][2] = 1; C[1][2] = C[2][1] = 0; C[1][3] = C[3][1]=5; C[2][3]=C[3][2] =1;
   //Eigenvalues found using Matlab
   double * eigenvalues_exact = new double[n];
   eigenvalues_exact[0] = -2.000; eigenvalues_exact[1] = 0.7574;eigenvalues_exact[2] = 3.0;eigenvalues_exact[3] = 9.2426;

   int k, l;
   double max_elem;
   max_element(C, &k, &l, &max_elem, n); // Calls on function that finds maximum value of matrix A
   double epsilon = 10e-100;
   double maxmax = max_elem*max_elem;
   int number_of_transformations = n*n*n;
   int counter = 0;
   while (maxmax >= epsilon && counter <= number_of_transformations){

   func(k, l, n, C, B, S);     // Calculates the updated matrix B and the new transformation matrix S
      for (int i = 0; i<n; i++){
         for (int j=0; j<n; j++){
            C[i][j] = B[i][j];
         }
      }

      max_element(C, &k, &l, &max_elem, n);

      maxmax = max_elem*max_elem;
      counter += 1;
   }
   double * eigenvalues_approx =  new double [n];
   for (int i =0; i<n; i++){
      eigenvalues_approx[i] = C[i][i];

   }
   double eps = 0.00001;
   int i = 0;
   int j = 0;
   while (i < n){
      if (j<n){
         if (fabs(eigenvalues_exact[i]-eigenvalues_approx[j]<=eps)){
            cout << "Sucess! Eigenvalue " << i << " within " << eps << " of correct eigenvalue" << endl;
         i += 1;
         }
         else {j += 1;}
      }
      else{
         cout << "Failure! " << "Eigenvalue " << i << " not within " << eps << " of correct eigenvalue" << endl;
         i = n;}

   }
}





int main(int argc, char* argv[]){
   string filename = argv[1];
   int n = atoi(argv[2]);
   double omega_r = atof(argv[4]);
   double * rho = new double[n+1];
   double * V = new double[n+1];

   string fileout = filename;


   // Defining dimensionless parameter rho
   double rho_0 = 0;
   double rho_n = atof(argv[3]);
   double h = (rho_n-rho_0)/(double)(n+1);  // Step length
   double hh = h*h;
   double const_diag = 2.0/hh;     // Constant term on diagonal
   double const_offdiag = -1.0/(hh);  // Constant non-diagonal elements
   for (int i=1; i<n+1; i++){rho[i] = rho_0 + double (i*h);}

   // Defining potential V
   for (int i = 1; i<n+1; i++){
      double rho_i = rho[i];
      V[i] = rho_i*rho_i;
   }


   // Creating empty matrices A, B and S
   double ** A = new double* [n];
   double ** B = new double* [n];
   double ** S = new double* [n];
   for (int i=0; i<n; i++){
      B[i] = new double [n];
      S[i] = new double [n];
      A[i] = new double [n];}
   unit_test(B, S);

   // Setting up inital matrix A
   for (int i = 0; i<n; i++){
      for (int j = 0; j<n; j++){
         if (i==j)
            //{A[i][j] = const_diag+V[i+1];}  //Non-interacting case
         {A[i][j]=const_diag+omega_r*omega_r*rho[i]*rho[i]+1/rho[i];}
         else if (fabs(i-j)==1){A[i][j] = const_offdiag;}      // Off-diagonal elements
         else {A[i][j] = 0;}
      }
   }

   // Setting up initial rotation matrix S
   for (int i = 0; i<n; i++){
      for (int j = 0; j<n; j++){
        if (i == j){S[i][j] = 1;}
        else {S[i][j] = 0;}
      }
   }




   // Check if matrix A was constructed correctly
   ofile.open(fileout.c_str());
   ofile << setiosflags(ios::showpoint | ios::uppercase);
   for (int i =0; i<n; i++){
      for (int j=0; j<n; j++){
         ofile << setw(12) << A[i][j];
      }
      ofile << endl;

   }
   ofile.close();



   int k, l;
   double max_elem;
   max_element(A, &k, &l, &max_elem, n); // Calls on function that finds maximum value of matrix A
   double epsilon = 10.0e-10;
   double maxmax = max_elem*max_elem;
   double number_of_transformations = (double)n*(double)n*(double)n;
   int counter = 0;
   while (maxmax >= epsilon && counter <= number_of_transformations){

      func(k, l, n, A, B, S);     // Calculates the updated matrix B and the new transformation matrix S
      for (int i = 0; i<n; i++){
         for (int j=0; j<n; j++){
            A[i][j] = B[i][j];
         }
      }

      max_element(A, &k, &l, &max_elem, n);

      maxmax = max_elem*max_elem;
      counter += 1;


   }
   for (int i =0; i<n; i++){
     for (int j=0; j<n; j++){
       if (i!= j){ B[i][j] = 0;}
     }
   }


   int one = 0;
   int two = 0;
   int three = 0;
   double l1 = 100;
   double l2 = 100;
   double l3 = 100;
   double l4 = 100;
   for (int i = 0; i<n; i++){
      for (int j = 0; j<n; j++){
         double bi = B[i][j];
         if (bi < l1 && bi > 0){
            l1 = bi;
            one = i;}
         else if (bi > l1 && bi < l2){
            l2 = bi;
            two = i;}
         else if (bi > l2 && bi < l3){
            l3 = bi;
            three = i;}
         else if (bi > l3 && bi < l4){l4 = bi;}
      }
   }
   /*
   ofile.open(fileout.c_str());
   ofile << setiosflags(ios::showpoint | ios::uppercase);
   for (int j=0; j<n; j++){
       double S1j = S[j][one]*S[j][one];
       double S2j = S[j][two]*S[j][two];
       double S3j = S[j][three]*S[j][three];
       ofile << setw(15) << S1j << setw(15) << S2j << setw(15) << S3j << endl;;
   }
   */


   ofile.close();

   cout << counter << endl;
   cout << l1 << endl;
   cout << l2 << endl;
   cout << l3 << endl;
   cout << l4 << endl;


   delete [] A; delete [] B; delete []S; delete []rho; delete [] V;
   return 0;
}
