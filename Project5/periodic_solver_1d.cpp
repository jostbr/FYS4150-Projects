
#include "periodic_solver_1d.hpp"

/* Constructor calling constructor of base class with input parameters. */
periodic_solver_1d::periodic_solver_1d(double dx, double dy,
                        int N, double T, std::string fileout) : rossby_solver_1d(dx, dy, N, T, fileout) {
    // No more initialization in derived class needed
}


/* Function that sets the initial psi and stores it as a member psi_0. */
void periodic_solver_1d::set_initial_condition(double* init_psi, double *init_zeta){
    /* Abort if users IC violates the periodic BC's. */
    if (init_psi[0] != init_psi[this->N-1]){
        std::cout << "Error: Chosen initial condition does not satisfy periodic BC!" << std::endl;
        std::cout << "Terminating program..." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->N; i++){
        this->psi_0[i] = init_psi[i];       // Store initial psi in member variable
        this->zeta_0[i] = init_zeta[i];     // Store initial zeta in member variable
    }
}


void periodic_solver_1d::periodic_euler(){
    // Euler time-stepping in periodic domain
}


void periodic_solver_1d::periodic_leapfrog(){
    // Leapfrog time-stepping in periodic domain
}


periodic_solver_1d::~periodic_solver_1d(){
    // No need for cleanup
}
