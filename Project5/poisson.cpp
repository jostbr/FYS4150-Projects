
# include "poisson.hpp"

/* Function implementing the Thomas algorithm for solving the linear system Ax = y where A
 * is a general tridiagonal matrix (NxN) with lower diag a (length N-1), main diag b (length N)
 * and upper diag c (length N-1). Further y is the known right-hand-side vector and solution
 * is the array to hold the solution x. Can e.g. be applied to solve the 1D Poisson equation
 *
 * d^2(solution)/dx^2 = y
 *
 * Note that this only solves for the interior points. */
void tridiag_general(double* a, double* b, double* c, double* y, int N, double*& solution){
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
void tridiag_ferrari(double* b, double* y, int N, double*& solution){
    //b[0] = 2.0;

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
void jacobi(double** g, double bc_0y, double bc_1y, double bc_x0, double bc_x1, double dx,
            double dy, int num_iter, int N_x, int N_y, double**& f){
    /* UNTESTED FUNCTION! */
    double dxdx = dx*dx;
    double dydy = dy*dy;
    double dxdxdydy = dxdx*dydy;
    double dxdx_pluss_dydy_2 = 2*(dxdx + dydy);

    for (int j = 0; j < N_y; j++){
        f[0][j] = bc_0y;         // Boundary condition at (x = 0, y)
        f[N_x-1][j] = bc_1y;     // Boundary condition at (x = 1, y)
    }

    for (int i = 0; i < N_y; i++){
        f[i][0] = bc_x0;         // Boundary condition at (x, y = 0)
        f[i][N_y-1] = bc_x1;     // Boundary condition at (x, y = 1)
    }

    for (int iter = 0; iter < num_iter; iter++){    // Iterations for convergence
        for (int i = 1; i < N_x-1; i++){
            for (int j = 1; j < N_y-1; j++){
                f[i][j] = (dydy*(f[i+1][j] + f[i-1][j]) + dxdx*(f[i][j+1] + f[i][j-1]) -
                        dxdxdydy*g[i][j])/dxdx_pluss_dydy_2;
            }
        }
    }
}
