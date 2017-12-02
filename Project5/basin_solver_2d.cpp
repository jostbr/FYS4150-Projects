
# include "basin_solver_2d.hpp"

basin_solver_2d::basin_solver_2d(double dx, double dy, int N_x, int N_y, double dt,
                     double T, std::string fileout) : rossby_solver(dx, dy, N_x, N_y, dt, T, fileout){
    alloc_array_1D(this->bc_0y, this->N_y);
    alloc_array_1D(this->bc_1y, this->N_y);
    alloc_array_1D(this->bc_x0, this->N_x);
    alloc_array_1D(this->bc_x1, this->N_x);

    for (int j = 0; j < this->N_y; j++){
        this->bc_0y[j] = 0.0;       // Zero flow through western boundary
        this->bc_1y[j] = 0.0;       // Zero flow through eastern boundary
    }

    for (int i = 0; i < this->N_x; i++){
        this->bc_x0[i] = 0.0;       // Zero flow through southern boundary
        this->bc_x1[i] = 0.0;       // Zero flow through northern boundary
    }
}


/* Function that sets the initial psi and zeta and stores it as members psi_0 and zeta_0. The
 * function also makes sure that the provided IC's satisfy the earlier supplied BC's. */
void basin_solver_2d::set_initial_condition(double* init_psi, double* init_zeta){
    double eps = 1.0E-10;   // Some small number to check for BC satisfaction

    /* Abort execution if users IC violates users BC's. */
    for (int j = 0; j < this->N_y; j++){
        if (fabs(init_psi[0 + j] - this->bc_0y[j]) > eps ||
                fabs(init_psi[(this->N_x-1)*this->N_y + j] - this->bc_1y[j]) > eps){
            std::cout << init_psi[0] << ", " << init_psi[this->N_x-1] << std::endl;
            std::cout << "Error: Chosen initial condition does not satisfy BC!" << std::endl;
            std::cout << "Terminating program..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < this->N_x; i++){
        if (fabs(init_psi[i*this->N_y + 0] - this->bc_x0[i]) > eps ||
                fabs(init_psi[i*this->N_y + this->N_y-1] - this->bc_x1[i]) > eps){
            std::cout << init_psi[0] << ", " << init_psi[this->N_x-1] << std::endl;
            std::cout << "Error: Chosen initial condition does not satisfy BC!" << std::endl;
            std::cout << "Terminating program..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < this->N_x; i++){
        for (int j = 0; j < this->N_y; j++){
            this->psi_0[i*this->N_y + j] = init_psi[i*this->N_y + j];       // Store initial psi in member variable
            this->zeta_0[i*this->N_y + j] = init_zeta[i*this->N_y + j];     // Store initial zeta in member variable
        }
    }

}


/* Function that solves the 2D vorticity equation in a bounded domain as a set of two coupled PDE's. In
 * particular, an advection equation is advanced forwrd in time and for every time-step, we solve the 2D
 * Poisson equation using an iterative solver (Jacobi method). */
void basin_solver_2d::basin_leapfrog(){
    this->display_model_config();

    /* Initialize variables and allocate memory. */
    /* ================================================================================== */
    double alpha = this->dt/(2.0*this->dx);     // Pre-calculate constant for Euler step
    double gamma = this->dt/this->dx;           // Pre-calculate constant for leapfrog steps
    double t = 0.0;                             // To keep track of time
    int array_size_2d = this->N_x*this->N_y;

    double *psi_prev, *psi_curr, *zeta_pp, *zeta_prev, *zeta_curr;    // Prev-prev, prev and current states
    alloc_array_1D(psi_prev, array_size_2d);
    alloc_array_1D(psi_curr, array_size_2d);
    alloc_array_1D(zeta_prev, array_size_2d);
    alloc_array_1D(zeta_pp, array_size_2d);
    alloc_array_1D(zeta_curr, array_size_2d);

    /* Set initial condition and apply boundary conditions. */
    /* ================================================================================== */
    for (int i = 0; i < this->N_x; i++){
        for (int j = 0; j < this->N_y; j++){
            psi_prev[i*this->N_y + j] = this->psi_0[i*this->N_y + j]; // Set previous psi to initial psi (takes care of BC as well)
            zeta_pp[i*this->N_y + j] = this->zeta_0[i*this->N_y + j]; // Set prev-prev zeta to initial zeta (takes care of BC as well)

            psi_curr[0 + j] = this->bc_0y[j];                           // Also apply BC for curr arrays
            psi_curr[(this->N_x-1)*this->N_y + j] = this->bc_1y[j];
            zeta_curr[0 + j] = zeta_pp[0 + j];                          // Also apply BC for curr arrays
            zeta_curr[(this->N_x-1)*this->N_y + j] = zeta_pp[(this->N_x-1)*this->N_y + j];
        }

        psi_curr[i*this->N_y + 0] = this->bc_x0[i];             // Also apply BC for curr arrays
        psi_curr[i*this->N_y + this->N_y-1] = this->bc_x1[i];
        zeta_curr[i*this->N_y + 0] = zeta_pp[i*this->N_y + 0];  // Also apply BC for curr arrays
        zeta_curr[i*this->N_y + this->N_y-1] = zeta_pp[i*this->N_y + this->N_y-1];
    }

    this->write_state_to_file(0.0, psi_prev);   // Write initial condition to file

    /* Initial Euler step. */
    /* ================================================================================== */
    for (int i = 1; i < this->N_x-1; i++){
        for (int j = 1; j < this->N_y-1; j++){
            zeta_prev[i*this->N_y + j] = this->zeta_0[i*this->N_y + j] -
                    alpha*(this->psi_0[(i+1)*this->N_y + j] - this->psi_0[(i-1)*this->N_y + j]);
        }
    }

    poisson_jacobi(zeta_prev, this->bc_0y, this->bc_1y, this->bc_x0, this->bc_x1,
                   this->dx, this->dy, this->N_x, this->N_y, 50, psi_prev);     // Solve 2D poisson eq. to get psi

    this->write_state_to_file(0.0, psi_prev);

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 2; t < T; n++){
        /* STEP 1: Advance vorticity forward in time (leapfrog). */
        for (int i = 1; i < this->N_x-1; i++){
            for (int j = 1; j < this->N_y-1; j++){
                zeta_curr[i*this->N_y + j] = zeta_pp[i*this->N_y + j] -
                        gamma*(psi_prev[(i+1)*this->N_y + j] - psi_prev[(i-1)*this->N_y + j]);
            }
        }

        /* STEP 2: Solve the 2D Poisson equation to update streamfunction. */
        poisson_jacobi(zeta_curr, this->bc_0y, this->bc_1y, this->bc_x0, this->bc_x1,
                       this->dx, this->dy, this->N_x, this->N_y, 50, psi_curr);

        for (int i = 1; i < this->N_x-1; i++){
            for (int j = 1; j < this->N_y-1; j++){
                psi_prev[i*this->N_y + j] = psi_curr[i*this->N_y + j];      // Set prev psi to curr for next time-step
                zeta_pp[i*this->N_y + j] = zeta_prev[i*this->N_y + j];      // Set prev-prev zeta to prev for next time-step
                zeta_prev[i*this->N_y + j] = zeta_curr[i*this->N_y + j];    // Set prev zeta to curr for next time-step
            }
        }

        t += this->dt;

        if (n % 200 == 0){
            this->write_state_to_file(t, psi_curr);
        }
    }

    free_array_1D(psi_prev);
    free_array_1D(psi_curr);
    free_array_1D(zeta_pp);
    free_array_1D(zeta_prev);
    free_array_1D(zeta_curr);
}


basin_solver_2d::~basin_solver_2d(){
    free_array_1D(this->bc_0y);
    free_array_1D(this->bc_1y);
    free_array_1D(this->bc_x0);
    free_array_1D(this->bc_x1);
}
