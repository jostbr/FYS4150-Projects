
#include "basin_solver_1d.hpp"

/* Constructor calling constructor of base class with input parameters. */
basin_solver_1d::basin_solver_1d(double dx, int N,
                     double dt, double T, std::string fileout) : rossby_solver(dx, N, dt, T, fileout) {
    this->bc_0 = 0.0;   // Basin boundary condition; psi(0,t) = 0
    this->bc_N = 0.0;   // Basin boundary condition; psi(1,t) = 0
}


/* Function that sets the initial psi and zeta and stores it as members psi_0 and zeta_0. The
 * function also makes sure that the provided IC's satisfy the earlier supplied BC's. */
void basin_solver_1d::set_initial_condition(double* init_psi, double* init_zeta){
    double eps = 1.0E-10;   // Some small numbner to check for BC satisfaction

    /* Abort execution if users IC violates users BC's. */
    if (fabs(init_psi[0] - this->bc_0) > eps || fabs(init_psi[this->N-1] - this->bc_N) > eps){
        std::cout << init_psi[0] << ", " << init_psi[this->N-1] << std::endl;
        std::cout << "Error: Chosen initial condition does not satisfy BC!" << std::endl;
        std::cout << "Terminating program..." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->N; i++){
        this->psi_0[i] = init_psi[i];       // Store initial psi in member variable
        this->zeta_0[i] = init_zeta[i];     // Store initial zeta in member variable
    }
}


/* Function that solves the vorticity equation as a set of two coupled PDE's. Uses Forward Euler
 * time-stepping for the vorticity and centered differences for the the streamfunction. Solves the
 * 1D Poisson equation at every time step using a tridiagonal linear algebra solver (Thomas algorithm)*/
void basin_solver_1d::basin_euler(){
    this->display_model_config();

    /* Initialize variables and allocate memory. */
    /* ================================================================================== */
    double alpha = this->dt/(2.0*this->dx);     // Pre-calculate constant for Euler step
    double dxdx = this->dx*this->dx;            // Pre-calculate for use in Poisson eq. below
    double t = 0.0;                             // To keep track of time

    double *psi_prev, *psi_curr, *zeta_prev, *zeta_curr;    // Previous and current states
    alloc_array_1D(psi_prev, this->N);
    alloc_array_1D(psi_curr, this->N);
    alloc_array_1D(zeta_prev, this->N);
    alloc_array_1D(zeta_curr, this->N);

    double *diag, *rhs_tridiag;                 // Arrays needed or call to tridiag solver
    alloc_array_1D(diag, this->N-2);
    alloc_array_1D(rhs_tridiag, this->N-2);

    /* Apply boundary consditions. */
    /* ================================================================================== */
    for (int i = 0; i < this->N; i++){
        psi_prev[i] = this->psi_0[i];       // Set previous psi to initial psi (takes care of BC as well)
        zeta_prev[i] = this->zeta_0[i];     // Set previous zeta to initial zeta (takes care of BC as well)
    }

    this->write_state_to_file(0.0, psi_prev);   // Write initial condition to file

    psi_curr[0] = this->bc_0; psi_curr[this->N-1] = this->bc_N;                 // Also apply BC here
    zeta_curr[0] = zeta_prev[0]; zeta_curr[this->N-1] = zeta_prev[this->N-1];   // Also apply BC here

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 1; t < T; n++){
        /* STEP 1: Advance vorticity for ward in time (forward Euler). */
        for (int i = 1; i < this->N-1; i++){
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1]);
        }

        /* STEP 2: Solve the 1D Poisson equation to update streamfunction. */
        for (int i = 1; i < this->N-1; i++){    // Loop over all interior zeta-values
            rhs_tridiag[i-1] = -dxdx*zeta_curr[i];     // Right-hand-side of Poisson eq.
        }

        tridiag_ferrari(diag, rhs_tridiag, this->N-2, psi_curr+1);   // Solve Poisson eq. in interior


        for (int i = 1; i < this->N-1; i++){
            psi_prev[i] = psi_curr[i];      // Set previous psi to current for next iter
            zeta_prev[i] = zeta_curr[i];    // Set previous zeta to current for next iter
        }

        t += this->dt;

        if (n % 200 == 0){
            this->write_state_to_file(t, psi_curr);
        }
    }

    free_array_1D(psi_prev);
    free_array_1D(psi_curr);
    free_array_1D(zeta_prev);
    free_array_1D(zeta_curr);
    free_array_1D(diag);
    free_array_1D(rhs_tridiag);
}

/* Function that solves the vorticity equation as a set of two coupled PDE's. Uses
 * Leapfrog time-stepping for the vorticity (forward Euler for the firat time step)
 * and centered differences for the the streamfunction. Solves the 1D Poisson equation
 * at every time step using a tridiagonal linear algebra solver (Thomas algorithm)*/
void basin_solver_1d::basin_leapfrog(){
    this->display_model_config();

    /* Initialize variables and allocate memory. */
    /* ================================================================================== */
    double alpha = this->dt/(2.0*this->dx);     // Pre-calculate constant for Euler step
    double gamma = this->dt/this->dx;           // Pre-calculate constant for leapfrog steps
    double dxdx = this->dx*this->dx;            // Pre-calculate for use in Poisson eq. below
    double t = 0.0;                             // To keep track of time

    double *psi_prev, *psi_curr, * zeta_prev, *zeta_pp, *zeta_curr;    // Prev-prev, prev and current states
    alloc_array_1D(psi_prev, this->N);
    alloc_array_1D(psi_curr, this->N);
    alloc_array_1D(zeta_prev, this->N);
    alloc_array_1D(zeta_pp, this->N);
    alloc_array_1D(zeta_curr, this->N);

    double *diag, *rhs_tridiag;                 // Arrays needed or call to tridiag solver
    alloc_array_1D(diag, this->N-2);
    alloc_array_1D(rhs_tridiag, this->N-2);

    /* Apply boundary consditions. */
    /* ================================================================================== */
    for (int i = 0; i < this->N; i++){
        psi_prev[i] = this->psi_0[i];       // Set previous psi to initial psi (takes care of BC as well)
        zeta_pp[i] = this->zeta_0[i];     // Set prev-prev zeta to initial zeta (takes care of BC as well)
    }

    this->write_state_to_file(0.0, psi_prev);   // Write initial condition to file

    psi_curr[0] = this->bc_0; psi_curr[this->N-1] = this->bc_N;             // Also apply BC for curr arrays
    zeta_curr[0] = zeta_pp[0]; zeta_curr[this->N-1] = zeta_pp[this->N-1];   // Also apply BC for curr arrays

    /* Initial Euler step. */
    /* ================================================================================== */
    for (int i = 1; i < this->N-1; i++){
        zeta_prev[i] = this->zeta_0[i] - alpha*(this->psi_0[i+1] - this->psi_0[i-1]);
    }

    for (int i = 1; i < this->N-1; i++){    // Loop over all interior zeta-values
        rhs_tridiag[i-1] = -dxdx*zeta_curr[i];     // Right-hand-side of Poisson eq.
    }

    tridiag_ferrari(diag, rhs_tridiag, this->N-2, psi_prev+1);   // Solve Poisson eq. in interior

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 2; t < T; n++){
        /* STEP 1: Advance vorticity for ward in time (leapfrog). */
        for (int i = 1; i < this->N-1; i++){
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1]);
        }

        /* STEP 2: Solve the 1D Poisson equation to update streamfunction. */
        for (int i = 1; i < this->N-1; i++){    // Loop over all interior zeta-values
            rhs_tridiag[i-1] = -this->dx*this->dx*zeta_curr[i];     // Right-hand-side of Poisson eq.
        }

        tridiag_ferrari(diag, rhs_tridiag, this->N-2, psi_curr+1);   // Solve Poisson eq. in interior


        for (int i = 1; i < this->N-1; i++){
            psi_prev[i] = psi_curr[i];      // Set previous psi to current for next iter
            zeta_pp[i] = zeta_prev[i];      // Set prev-prev zeta to prev for next iter
            zeta_prev[i] = zeta_curr[i];    // Set previous zeta to current for next iter
        }

        t += this->dt;

        if (n % 200 == 0){
            this->write_state_to_file(t, psi_curr);
        }
    }

    free_array_1D(psi_prev);
    free_array_1D(psi_curr);
    free_array_1D(zeta_prev);
    free_array_1D(zeta_pp);
    free_array_1D(zeta_curr);
    free_array_1D(diag);
    free_array_1D(rhs_tridiag);
}


basin_solver_1d::~basin_solver_1d(){
    // No need for cleanup
}
