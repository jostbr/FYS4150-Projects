
#ifndef POISSON_HPP
#define POISSON_HPP

# include <iostream>
# include <cmath>
# include "array_alloc.hpp"

void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution);
void tridiag_ferrari(double* b, double* y, int N, double* solution);
void poisson_jacobi(double* g, double* bc_0y, double* bc_1y, double* bc_x0, double* bc_x1, double dx,
            double dy, int N_x, int N_y, int max_iter, double* f);

#endif // POISSON_HPP
