#ifndef POISSON_HPP
#define POISSON_HPP

void tridiag(double* a, double* b, double* c, double* y, int N, double* solution);
void jacobi(double** g, double bc_0y, double bc_1y, double bc_x0, double bc_x1, double dx,
            double dy, int num_iter, int N_x, int N_y, double**& f);

#endif // POISSON_HPP
