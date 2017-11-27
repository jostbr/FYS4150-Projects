
#include "periodic_solver_1d.hpp"

/* Constructor calling constructor of base class with input parameters. */
periodic_solver_1d::periodic_solver_1d(double dx, int N_x,
                        double dt, double T, std::string fileout) : rossby_solver(dx, N_x, dt, T, fileout) {
    // No more initialization in derived class needed
}


/* Function that sets the initial psi and stores it as a member psi_0. */
void periodic_solver_1d::set_initial_condition(double* init_psi, double *init_zeta){
    double eps = 1.0E-10;   // Some small numbner to check for BC satisfaction
    /* Abort if users IC violates the periodic BC's. */
    if ((fabs(init_psi[0] - init_psi[this->N_x-1]) > eps) || (fabs(init_zeta[0] - init_zeta[this->N_x-1]) > eps)){
        std::cout << "Error: Chosen initial condition does not satisfy periodic BC!" << std::endl;
        std::cout << "Terminating program..." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->N_x; i++){
        this->psi_0[i] = init_psi[i];       // Store initial psi in member variable
        this->zeta_0[i] = init_zeta[i];     // Store initial zeta in member variable
    }
}


/* Function that solves the vorticity equation with periodic boundary consitions as a set of two coupled PDE's.
 * Uses Forward Euler time-stepping for the vorticity and centered differences for the the streamfunction. Solves
 * the 1D Poisson equation at every time step using a Armadillo LU decomposition linear algebra solver */
void periodic_solver_1d::periodic_euler(){
    this->display_model_config();

    /* Initialize variables and allocate memory. */
    /* ================================================================================== */
    double alpha = this->dt/(2.0*this->dx);     // Pre-calculate constant for Euler step
    double dxdx = this->dx*this->dx;            // Pre-calculate for use in Poisson eq. below
    double t = 0.0;                             // To keep track of time

    arma::vec psi_prev(this->N_x);
    arma::vec psi_curr(this->N_x);
    arma::vec zeta_prev(this->N_x);
    arma::vec zeta_curr(this->N_x);

    arma::vec rhs_poisson(this->N_x);     // Vector to hold right-hand-side of Poisson eq.
    arma::mat A = arma::zeros(this->N_x, this->N_x);    // Matrix to store second-derivative coefficients

    this->initialize_periodic_matrix(A);    // Initialize A for Armadillo LU solver used every time-step

    /* Set initial conditions. */
    /* ================================================================================== */
    for (int i = 0; i < this->N_x; i++){
        psi_prev[i] = this->psi_0[i];       // Set previous psi to initial psi (takes care of BC as well)
        zeta_prev[i] = this->zeta_0[i];     // Set previous zeta to initial zeta (takes care of BC as well)
    }

    this->write_state_to_file(0.0, psi_prev.memptr());   // Write initial condition to file

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 1; t < T; n++){
        /* STEP 1: Advance vorticity forward in time (forward Euler). */
        for (int i = 1; i < this->N_x-1; i++){
            zeta_curr[i] = zeta_prev[i] - alpha*(psi_prev[i+1] - psi_prev[i-1]);
        }

//        zeta_curr[0] = zeta_prev[0] - alpha*(psi_prev[1] - psi_prev[this->N_x-1]);                    // Periodic BC
//        zeta_curr[this->N_x-1] = zeta_prev[this->N_x-1] - alpha*(psi_prev[0] - psi_prev[this->N_x-2]);    // Periodic BC

        zeta_curr[0] = zeta_prev[0] - alpha*(psi_prev[1] - psi_prev[this->N_x-2]);                    // Periodic BC
        zeta_curr[this->N_x-1] = zeta_prev[this->N_x-1] - alpha*(psi_prev[1] - psi_prev[this->N_x-2]);    // Periodic BC

        /* STEP 2: Solve the 1D Poisson equation to update streamfunction. */
        for (int i = 0; i < this->N_x; i++){    // Loop over all interior zeta-values
            rhs_poisson[i] = -dxdx*zeta_curr[i];     // Right-hand-side of Poisson eq. (interior + boundary)
        }

        psi_curr = arma::solve(A, rhs_poisson);     // Solve 1D Poisson eq. to get psi at current time-step

        for (int i = 0; i < this->N_x; i++){
            psi_prev[i] = psi_curr[i];      // Set previous psi to current for next iter
            zeta_prev[i] = zeta_curr[i];    // Set previous zeta to current for next iter
        }

        t += this->dt;

        if (n % 200 == 0){
            this->write_state_to_file(t, psi_curr.memptr());
        }
    }
}


/* Function that solves the vorticity equation with periodic boundary conditions as a set of two
 * coupled PDE's. Uses Leapfrog time-stepping for the vorticity (forward Euler for the firat time
 * step) and centered differences for the the streamfunction. Solves the 1D Poisson equation
 * at every time step using a Armadillo LU decomposition linear algebra solver. */
void periodic_solver_1d::periodic_leapfrog(){
    this->display_model_config();

    /* Initialize variables and allocate memory. */
    /* ================================================================================== */
    double alpha = this->dt/(2.0*this->dx);     // Pre-calculate constant for Euler step
    double gamma = this->dt/this->dx;           // Pre-calculate constant for leapfrog steps
    double dxdx = this->dx*this->dx;            // Pre-calculate for use in Poisson eq. below
    double t = 0.0;                             // To keep track of time

    arma::vec psi_prev(this->N_x);
    arma::vec psi_curr(this->N_x);
    arma::vec zeta_pp(this->N_x);
    arma::vec zeta_prev(this->N_x);
    arma::vec zeta_curr(this->N_x);

    arma::vec rhs_poisson(this->N_x);     // Vector to hold right-hand-side of Poisson eq.
    arma::mat A = arma::zeros(this->N_x, this->N_x);    // Matrix to store second-derivative coefficients

    this->initialize_periodic_matrix(A);    // Initialize A for Armadillo LU solver used every time-step

    /* Set initial conditions. */
    /* ================================================================================== */
    for (int i = 0; i < this->N_x; i++){
        psi_prev[i] = this->psi_0[i];       // Set previous psi to initial psi (takes care of BC as well)
        zeta_pp[i] = this->zeta_0[i];     // Set previous zeta to initial zeta (takes care of BC as well)
    }

    this->write_state_to_file(0.0, psi_prev.memptr());   // Write initial condition to file

    /* Initial Euler step. */
    /* ================================================================================== */
    for (int i = 1; i < this->N_x-1; i++){
        zeta_prev[i] = this->zeta_0[i] - alpha*(this->psi_0[i+1] - this->psi_0[i-1]);
    }

    zeta_prev[0] = this->zeta_0[0] - alpha*(this->psi_0[1] - this->psi_0[this->N_x-2]);                    // Periodic BC
    zeta_prev[this->N_x-1] = this->zeta_0[this->N_x-1] - alpha*(this->psi_0[1] - this->psi_0[this->N_x-2]);    // Periodic BC

    for (int i = 0; i < this->N_x; i++){    // Loop over all interior zeta-values
        rhs_poisson[i] = -dxdx*zeta_prev[i];     // Right-hand-side of Poisson eq. (interior + boundary)
    }

    psi_prev = arma::solve(A, rhs_poisson);     // Solve 1D Poisson eq. to get psi at time-step 1

    /* Main loop over time using leapfrog time-stepping for vorticity. */
    /* ================================================================================== */
    for (int n = 2; t < T; n++){
        /* STEP 1: Advance vorticity for ward in time (leapfrog). */
        for (int i = 1; i < this->N_x-1; i++){
            zeta_curr[i] = zeta_pp[i] - gamma*(psi_prev[i+1] - psi_prev[i-1]);
        }

        zeta_curr[0] = zeta_pp[0] - gamma*(psi_prev[1] - psi_prev[this->N_x-2]);                    // Periodic BC
        zeta_curr[this->N_x-1] = zeta_pp[this->N_x-1] - gamma*(psi_prev[1] - psi_prev[this->N_x-2]);    // Periodic BC

        /* STEP 2: Solve the 1D Poisson equation to update streamfunction. */
        for (int i = 0; i < this->N_x; i++){    // Loop over all interior zeta-values
            rhs_poisson[i] = -dxdx*zeta_curr[i];     // Right-hand-side of Poisson eq. (interior + boundary)
        }

        psi_curr = arma::solve(A, rhs_poisson);     // Solve 1D Poisson eq. to get psi at current time-step


        for (int i = 0; i < this->N_x; i++){
            psi_prev[i] = psi_curr[i];      // Set previous psi to current for next iter
            zeta_pp[i] = zeta_prev[i];      // Set prev-prev zeta to prev for next iter
            zeta_prev[i] = zeta_curr[i];    // Set previous zeta to current for next iter
        }

        t += this->dt;

        if (n % 200 == 0){
            this->write_state_to_file(t, psi_curr.memptr());
        }
    }
}


/* Function to initialize periodic BC matrix to pass to Armadillo to solve Poisson at every time step. */
void periodic_solver_1d::initialize_periodic_matrix(arma::mat& A){
    for (int i = 0; i < this->N_x; i++){
        for (int j = 0; j < this->N_x; j++){
            if (i == j){
                A(i,j) = 2.0;       // Main diagonal elements
            }

            else if (fabs(i - j) == 1){
                A(i,j) = -1.0;      // Lower and upper diagonal elements
            }
        }
    }

    A(0,this->N_x-2) = -1.0;      // Semi-corners value due to periodic BC
    A(this->N_x-1,1) = -1.0;      // Semi-corners value due to periodic BC
}


periodic_solver_1d::~periodic_solver_1d(){
    // No need for cleanup
}
