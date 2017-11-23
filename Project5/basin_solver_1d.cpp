
#include "basin_solver_1d.hpp"

/* Constructor calling constructor of base class wit input parameters. */
basin_solver_1d::basin_solver_1d(double dx, double dt,
                     int N, double T) : rossby_solver_1d(dx, dt, N, T) {
    // No more initialization in derived class needed
}


/* Function that sets what psi should be at the boundaries. */
void basin_solver_1d::set_boundary_conditions(double bc_0, double bc_N){
    this->bc_0 = bc_0;
    this->bc_N = bc_N;
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
    }

    for (int i = 0; i < this->N; i++){
        this->psi_0[i] = init_psi[i];
        this->zeta_0[i] = init_zeta[i];
    }
}


/* Function that solves the vorticity equation as a coupled set of two PDE's. Uses Forward Euler
 * time-stepping for the vorticity and centered differences for the the streamfunction. */
void basin_solver_1d::basin_euler(){
    double alpha = this->dt/(2.0*this->dx);     // Pre-calcualte constant for scheme below
    double t = 0.0;                             // To keep track of time

    double *psi_prev, *psi_curr, *zeta_prev, *zeta_curr;    // Prevois and current states
    alloc_array_1D(psi_prev, this->N);
    alloc_array_1D(psi_curr, this->N);
    alloc_array_1D(zeta_prev, this->N);
    alloc_array_1D(zeta_curr, this->N);

    double *diag, *rhs_tridiag;                 // Arrays needed or call to tridiag solver
    alloc_array_1D(diag, this->N-2);
    alloc_array_1D(rhs_tridiag, this->N-2);

    for (int i = 0; i < this->N; i++){
        psi_prev[i] = this->psi_0[i];       // Set previous psi to initial psi (takes care of BC as well)
        zeta_prev[i] = this->zeta_0[i];     // Set previous zeta to initial zeta (takes care of BC as well)
    }

    psi_curr[0] = this->bc_0; psi_curr[this->N-1] = this->bc_N;                       // Also apply BC here
    zeta_curr[0] = zeta_prev[0]; zeta_curr[this->N-1] = zeta_prev[this->N-1];   // Also apply BC here

    while (t < this->T){

        /* STEP 1: Time-step for vorticity. */
        for (int i = 1; i < this->N-1; i++){
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1]);
        }

        /* STEP 2: Solve the 2D Poisson equation. */
        for (int i = 1; i < this->N-1; i++){
            rhs_tridiag[i-1] = -this->dx*this->dx*zeta_curr[i];
        }

        tridiag_ferrari(diag, rhs_tridiag, this->N-2, psi_curr+1);   // Solve Poisson eq. in interior


        for (int i = 1; i < this->N-1; i++){    // All points except at boundaries (always)
            psi_prev[i] = psi_curr[i];      // Set previous psi to current for next iter
            zeta_prev[i] = zeta_curr[i];    // Set previous zeta to current for next iter
        }

        t += this->dt;
        // Write state to file every some t step.
    }

    free_array_1D(psi_prev);
    free_array_1D(psi_curr);
    free_array_1D(zeta_prev);
    free_array_1D(zeta_curr);
    free_array_1D(diag);
    free_array_1D(rhs_tridiag);
    std::cout << "Sup?" << std::endl;
}


void basin_solver_1d::basin_leapfrog(){
    // Leapfrog time-stepping in basin domain
}


basin_solver_1d::~basin_solver_1d(){
    // No need for cleanup
}
