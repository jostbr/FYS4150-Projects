#ifndef JACOBI_H
#define JACOBI_H

void jacobi_master(arma::mat A, int N);
double get_max_non_diag(arma::mat A, int N, int* k, int* l);

#endif // JACOBI_H
