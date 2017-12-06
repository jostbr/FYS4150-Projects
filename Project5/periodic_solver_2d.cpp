
# include "periodic_solver_2d.hpp"

periodic_solver_2d::periodic_solver_2d(double dx, double dy, int N_x, int N_y, double dt,
                     double T, std::string fileout) : rossby_solver(dx, dy, N_x, N_y, dt, T, fileout){
    // No need for more initialization
}


/* Function that sets the initial psi and zeta and stores it as members psi_0 and zeta_0. The
 * function also makes sure that the provided IC's satisfy the earlier supplied BC's. */
void periodic_solver_2d::set_initial_condition(double* init_psi, double* init_zeta){
    double eps = 1.0E-10;   // Some small number to check for BC satisfaction

    /* Abort execution if users IC violates users BC's. */
    for (int j = 0; j < this->N_y; j++){
        if (fabs(init_psi[0 + j] - init_psi[(this->N_x-1)*this->N_y + j]) > eps){
            std::cout << init_psi[0] << ", " << init_psi[this->N_x-1] << std::endl;
            std::cout << "Error: Chosen initial condition does not satisfy BC!" << std::endl;
            std::cout << "Terminating program..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for (int i = 0; i < this->N_x; i++){
        if (fabs(init_psi[i*this->N_y + 0] - init_psi[i*this->N_y + (this->N_y-1)]) > eps){
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


/* Function that solves the 2D vorticity equation in a periodic domain as a set of two coupled PDE's. In
 * particular, an advection equation is advanced forwrd in time and for every time-step, we solve the 2D
 * Poisson equation using an iterative solver (Jacobi method). */
void periodic_solver_2d::periodic_leapfrog(){
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
        }
    }

    this->write_state_to_file(0.0, psi_prev);   // Write initial condition to file

    /* Initial Euler step. */
    /* ================================================================================== */
    for (int i = 1; i < this->N_x-1; i++){
        for (int j = 0; j < this->N_y; j++){
            zeta_prev[i*this->N_y + j] = this->zeta_0[i*this->N_y + j] -
                    alpha*(this->psi_0[(i+1)*this->N_y + j] - this->psi_0[(i-1)*this->N_y + j]);
        }
    }

    for (int j = 0; j < this->N_y; j++){
        zeta_prev[0*this->N_y + j] = this->zeta_0[0*this->N_y + j] -
                alpha*(this->psi_0[1*this->N_y + j] - this->psi_0[(this->N_x-2)*this->N_y + j]);    // Boundary x = 0
        zeta_prev[(this->N_x-1)*this->N_y + j] = zeta_prev[0*this->N_y + j];                        // Boundary x = 1
    }

    poisson_jacobi_periodic(zeta_prev, this->dx, this->dy, this->N_x, this->N_y, 50, psi_prev);     // Solve 2D poisson eq. to get psi

    this->write_state_to_file(0.0, psi_prev);

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 2; t < T; n++){
        /* STEP 1: Advance vorticity forward in time (leapfrog). */
        for (int i = 1; i < this->N_x-1; i++){
            for (int j = 0; j < this->N_y; j++){
                zeta_curr[i*this->N_y + j] = zeta_pp[i*this->N_y + j] -
                        gamma*(psi_prev[(i+1)*this->N_y + j] - psi_prev[(i-1)*this->N_y + j]);
            }
        }

        for (int j = 0; j < this->N_y; j++){
            zeta_curr[0*this->N_y + j] = zeta_pp[0*this->N_y + j] -
                    gamma*(psi_prev[1*this->N_y + j] - psi_prev[(this->N_x-2)*this->N_y + j]);  // Boundary x = 0
            zeta_curr[(this->N_x-1)*this->N_y + j] = zeta_curr[0*this->N_y + j];                // Boundary x = 1
        }

        /* STEP 2: Solve the 2D Poisson equation to update streamfunction. */
        poisson_jacobi_periodic(zeta_curr, this->dx, this->dy, this->N_x, this->N_y, 50, psi_curr);

        for (int i = 0; i < this->N_x; i++){
            for (int j = 0; j < this->N_y; j++){
                psi_prev[i*this->N_y + j] = psi_curr[i*this->N_y + j];      // Set prev psi to curr for next time-step
                zeta_pp[i*this->N_y + j] = zeta_prev[i*this->N_y + j];      // Set prev-prev zeta to prev for next time-step
                zeta_prev[i*this->N_y + j] = zeta_curr[i*this->N_y + j];    // Set prev zeta to curr for next time-step
            }
        }

        t += this->dt;

        if (n % 50 == 0){
            this->write_state_to_file(t, psi_curr);
        }
    }

    free_array_1D(psi_prev);
    free_array_1D(psi_curr);
    free_array_1D(zeta_pp);
    free_array_1D(zeta_prev);
    free_array_1D(zeta_curr);
}


periodic_solver_2d::~periodic_solver_2d(){
    // No cleanup needed
}
