
# include "unit_tests.hpp"

/* Unit test for the general tridiag solver. The test is based on using small arbitrary
 * tridiagonal matrix (3x3) with a arbitrary known right-hand- side. We compute the
 * solution by hand and compare with the output of tridiag(). */
void TEST_tridiag_general(){
    int N = 3;
    double *a, *b, *c, *y, *solution;
    alloc_array_1D(a, N-1);        // Lower diagonal
    alloc_array_1D(b, N);          // Main diagonal
    alloc_array_1D(c, N-1);    // Upper diagonal
    alloc_array_1D(y, N);          // Right hand side of Av = g
    alloc_array_1D(solution, N);   // To hold solution from algorithm

    /* Testing with an abitrary matrix for the general algorithm. */
    a[0] = a[1] = c[0] = c[1] = 1.0;
    b[0] = 1.0;
    b[1] = 2.0;
    b[2] = 3.0;
    y[0] = y[1] = y[2] = -1.0;      // Right hand side fortest of general

    tridiag_general(a, b, c, y, N, solution);   // Call the Thomas algorithm

    /* Examine results and see if results are as computed by hand. */
    if (solution[0] == -1.5 && solution[1] == 0.5 && solution[2] == -0.5){
        std::cout << "TEST PASSED: General tridiag solver" << std::endl;
    }

    else {
        std::cout << "TEST FAILED: General tridiag solver!" << std::endl;
    }

    free_array_1D(a);
    free_array_1D(b);
    free_array_1D(c);
    free_array_1D(y);
    free_array_1D(solution);
}


/* Unit test for the general tridiag solver. The test is based on using small arbitrary
 * tridiagonal matrix (3x3) with a arbitrary known right-hand- side. We compute the
 * solution by hand and compare with the output of tridiag(). */
void TEST_tridiag_ferrari(){
    int N = 3;
    double *diag, *rhs, *solution;
    alloc_array_1D(diag, N);          // Main diagonal
    alloc_array_1D(rhs, N);          // Right hand side of Av = g
    alloc_array_1D(solution, N);   // To hold solution from algorithm

    /* Testing with an abitrary matrix for the general algorithm. */
    rhs[0] = rhs[1] = rhs[2] = -1.0;      // Right hand side fortest of general

    tridiag_ferrari(diag, rhs, N, solution);   // Call the Thomas algorithm

    /* Examine results and see if results are as computed by hand. */
    if (solution[0] == -1.5 && solution[1] == 0.5 && solution[2] == -0.5){
        std::cout << "TEST PASSED: Ferrari tridiag solver" << std::endl;
    }

    else {
        std::cout << "TEST FAILED: Ferrari tridiag solver" << std::endl;
    }

    free_array_1D(diag);
    free_array_1D(rhs);
    free_array_1D(solution);
}

/* Unit test for testing if the basin initial condition are set correctly. */
void TEST_set_initial_condition_basin(){
    int N = 10;
    double bc = 12.0;
    double *init_psi, *init_zeta, *set_init_psi, *set_init_zeta;
    std::string fileout = "filename.txt";

    alloc_array_1D(init_psi, N);
    alloc_array_1D(init_zeta, N);
    alloc_array_1D(set_init_psi, N);
    alloc_array_1D(set_init_zeta, N);

    for (int i = 0; i < N; i++){
        init_psi[i] = bc;
        init_zeta[i] = 0.0;
    }

    basin_solver_1d solver(1.0, 1.0, N, 1.0, fileout);
    solver.set_boundary_conditions(bc, bc);
    solver.set_initial_condition(init_psi, init_zeta);
    solver.get_initial_condition(set_init_psi, set_init_zeta);

    for (int i = 0; i < N; i++){
        if (set_init_psi[i] != init_psi[i] || set_init_zeta[i] != init_zeta[i]){
            std::cout << "TEST FAILED: Basin initial condition" << std::endl;
            return;
        }
    }

    std::cout << "TEST PASSED: Basin initial condition" << std::endl;
}
