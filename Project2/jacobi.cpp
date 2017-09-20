#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;
ofstream ofile;

// Function to find the maximum off-diagonal element of a matrix A and the indices of this element.
void max_element(double ** A, int* k, int* l, double* max_elem, int n){
   *k = *l = *max_elem =0;     // Sets initial values to zero
    for (int i=0; i<n+1; i++){
       for (int j=i+1; j<n+1; j++){
	  double a = A[i][j];
          double aa = a*a;
          if (aa > *max_elem){
             *max_elem = A[i][j];
             *k = i;
             *l = j;
           }
        }
     }
}   

void func(double max_elem, int k, int l, int n, double ** A, double *** B, double ***S){
  
   double tau = (A[l][l]-A[k][k])/(2.0*A[k][l]);
   double t;
   double t_plus = -tau+sqrt(1+tau*tau);
   double t_minus = -tau-sqrt(1+tau*tau);
   // Finding smallest t
   if (t_plus < t_minus){t = t_plus;}
   else {t=t_minus;}
   double c = 1.0/sqrt(1+t*t);
   double s = t*c;
   double cc = c*c;
   double ss = s*s;
   double cs = c*s;
   cout << c << endl;

   for (int i=0; i<n+1; i++){
      for (int j=0; j<n+1; j++){
         if (i == j && i !=k && i != l){
           (*B)[i][j] = A[i][j];}
         else if (i!=k && i!=l && j==k){
            (*B)[i][j] = (*B)[j][i]= A[i][k]*c-A[i][l]*s;}
         else if (i!=k && i!=l && j == l){
            (*B)[i][j] = (*B)[j][i]= A[i][l]*c+A[i][k]*s;}
         else if (i==k && j==k){
            (*B)[i][j] = A[k][k]*cc-2.0*A[k][l]*cs + A[l][l]*ss;}
         else if (i==l && j == l){
            (*B)[i][j] = A[l][l]*cc+2*A[k][l]*cs+A[k][k]*ss;}
         else if (i==k && j==l){
            (*B)[i][j] = (*B)[j][i] = 0.0;}
         else {(*B)[i][j]=0.0;}
      }
   }
}   
   

int main(int argc, char* argv[]){
   string filename = argv[1];
   int n = atoi(argv[2]);
   double * rho = new double[n+1];
   double * V = new double[n+1];

   string fileout = filename;

   // Defining dimensionless paramter rho
   double rho_0 = 0; 
   double rho_n = 10e6;
   double h = (rho_n-rho_0)/(double)n;  // Step length
   double hh = h*h;
   double const_diag = 2.0/hh;     // Constant term on diagonal
   double const_offdiag = -1.0/(h*h);  // Constant non-diagonal elements
   rho[0] = rho_0;
   rho[n] = rho_n;
   for (int i=1; i<n; i++){rho[i] = rho_0 + i*h;}

   // Defining potential V
   for (int i = 0; i<n+1; i++){
      double rho_i = rho[i];
      V[i] = rho_i;
   }
      
   
   // Creating empty matrices A, B and S
   double ** A = new double* [n+1];
   double ** B = new double* [n+1];
   double ** S = new double* [n+1];
   for (int i=0; i<n+1; i++){
      B[i] = new double [n+1];
      S[i] = new double [n+1];
      A[i] = new double [n+1];}
   
   // Setting up inital matrix A
   for (int i = 0; i<n+1; i++){
      for (int j = 0; j<n+1; j++){
         if (i==j){A[i][j] = const_diag+V[i];}  //Diagonal elements
         else if (fabs(i-j)==1){A[i][j] = const_offdiag;}      // Off-diagonal elements
         else {A[i][j] = 0;}
      }
   }

    /* Check if matrix A was constructed correctly
   ofile.open(fileout.c_str());
   ofile << setiosflags(ios::showpoint | ios::uppercase);
   for (int i =0; i<n+1; i++){
      for (int j=0; j<n+1; j++){
         ofile << setw(12) << A[i][j];
      }
      ofile << endl;
      
   }
   ofile.close();
   */

   int k, l;
   double max;
   max_element(A, &k, &l, &max, n); // Calls on function that finds maximum value of matrix A
   double epsilon = 10e-40;
   double maxmax = max*max;
   int number_of_transformations = 10;
   int counter = 0;
   while (maxmax >= epsilon && counter <= number_of_transformations){
   
      func(max, k, l, n, A, &B, &S);     // Calculates the updated matrix B and the new transformation matrix S
      for (int i = 0; i<n+1; i++){
         for (int j=0; j<n+1; j++){
            A[i][j] = B[i][j];
         }
      }
      max_element(A, &k, &l, &max, n);
      maxmax = max*max;
      counter += 1;
      cout << k << "   " << l << endl;
      cout << max <<endl;
   }
   ofile.open(fileout.c_str());
   ofile << setiosflags(ios::showpoint | ios::uppercase);
   for (int i =0; i<n+1; i++){
      for (int j=0; j<n+1; j++){
         ofile << setw(12) << B[i][j];
      }
      ofile << endl;
      
   }
   ofile.close();
   cout << B[0][0] << "   "<< B[0][1]<< "   " << B[1][0]<< "    " << B[1][1] << endl;
   
   
  
   delete [] A; delete [] B; delete []S; delete []rho; delete [] V;
   return 0;
}
