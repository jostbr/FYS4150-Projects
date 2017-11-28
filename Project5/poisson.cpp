
# include "poisson.hpp"

/* Function implementing the Thomas algorithm for solving the linear system Ax = y where A
 * is a general tridiagonal matrix (NxN) with lower diag a (length N-1), main diag b (length N)
 * and upper diag c (length N-1). Further y is the known right-hand-side vector and solution
 * is the array to hold the solution x. Can e.g. be applied to solve the 1D Poisson equation
 *
 * d^2(solution)/dx^2 = y
 *
 * Note that this only solves for the interior points. */
void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution){
    /* STEP 1: Forward substitution. */
    for (int i = 1; i < N; i++){
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1];     // Eliminating lower diagonal
        y[i] = y[i] - (a[i-1]/b[i-1])*y[i-1];   // Corresponding change to RHS of eq.
    }

    /* STEP 2: Backward substitution. */
    solution[N-1] = y[N-1]/b[N-1];  // Special case for obtaining final element of solution

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] - c[i]*solution[i+1])/b[i];  // Eliminating upper diag and dividing by main diag
    }
}


/* Function for solving the linear system Ax = y where A is a special tridiagonal matrix (NxN) with
 * all elements along lower and upper diag equal to -1, while the main diag has all values equal to
 * 2. Further y is the known right-hand-side vector and solution is the array to hold the solution x.
 * Can e.g. be applied to solve the 1D Poisson equation
 *
 * d^2(solution)/dx^2 = y
 *
 * Note that this only solves for the interior points. */
void tridiag_ferrari(double* b, double* y, int N, double* solution){
    b[0] = 2.0;

    /* STEP 1: Forward substitution. */
    for (int i = 1; i < N; i++){
        b[i] = (i + 2)/((double)(i + 1));       // Eliminating lower diagonal
        y[i] = y[i] + (y[i-1]/b[i-1]);      // Corresponding change to RHS of eq.
    }

    /* STEP 2: Backward substitution. */
    solution[N-1] = y[N-1]/b[N-1];          // Special case for obtaining final element of solution

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] + solution[i+1])/b[i];  // Eliminating upper diag and dividing by main diag
    }
}


/* Function that implements the iterative Jacobi algorithm for solving the 2D Poisson equation
 *
 * d^2f/dx^2 + d^f/dy^2 = g
 *
 * with source function g(x,y) and (constant) boundary conditions
 *
 * f(0,y) = bc_0y
 * f(1,y) = bc_1y
 * f(x,0) = bc_x0
 * f(x,1) = bc_x1
 *
 * using num_iter number of iterations for convergence at each of the (N_x, N_y) spatial points. */
void poisson_jacobi(double* g, double* bc_0y, double* bc_1y, double* bc_x0, double* bc_x1, double dx,
            double dy, int N_x, int N_y, int max_iter, double* f){
    /* UNTESTED FUNCTION! */
    double dxdx = dx*dx;
    double dydy = dy*dy;
    double dxdxdydy = dxdx*dydy;
    double dxdx_pluss_dydy_2 = 2*(dxdx + dydy);

    for (int j = 0; j < N_y; j++){
        f[0 + j] = bc_0y[j];                // Boundary condition at (x = 0, y)
        f[(N_x-1)*N_y + j] = bc_1y[j];      // Boundary condition at (x = 1, y)
    }

    for (int i = 0; i < N_x; i++){
        f[i*N_y + 0] = bc_x0[i];            // Boundary condition at (x, y = 0)
        f[i*N_y + (N_y-1)] = bc_x1[i];      // Boundary condition at (x, y = 1)
    }

    int iter = 0;
    double diff = 1.0E+20;      // To check for convergence
    double eps = 1.0E-6;    // Tolerance for convergence
    double* f_tmp;          // To temporary hold solution for each iteration
    alloc_array_1D(f_tmp, N_x*N_y);

    /* Iterate until satisafactory convergence is reached (or max iter is reached). */
    while (iter <= max_iter && fabs(diff) > eps){
        diff = 0.0;

        for (int i = 0; i < N_x; i++){
            for (int j = 0; j < N_y; j++){
                f_tmp[i*N_y + j] = f[i*N_y + j];    // Need previous "solution"
            }
        }

        /* Do one sweep over the array step each point closer to the solution. */
        for (int i = 1; i < N_x-1; i++){
            for (int j = 1; j < N_y-1; j++){
                f[i*N_y + j] = (dydy*(f_tmp[(i+1)*N_y + j] + f_tmp[(i-1)*N_y + j]) +
                        dxdx*(f_tmp[i*N_y + (j+1)] + f_tmp[i*N_y + (j-1)]) -
                        dxdxdydy*g[i*N_y + j])/dxdx_pluss_dydy_2;
                diff += f[i*N_y + j] - f_tmp[i*N_y + j];
            }
        }

        //std::cout << diff << std::endl;

        iter++;
    }

    if (fabs(diff) > eps){
        //std::cout << "Did not reach satisfactory convergence!" << std::endl;
    }

    free_array_1D(f_tmp);
}
