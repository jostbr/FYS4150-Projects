
#include "periodic_solver_1d.hpp"

periodic_solver_1d::periodic_solver_1d(double dx, double dy,
                        int N, double T, std::string fileout) : rossby_solver_1d(dx, dy, N, T, fileout) {
    // No more initialization in derived class needed
}


/* Function that sets the initial psi and stores it as a member psi_0. */
void periodic_solver_1d::set_initial_condition(double* init_con){
    /* Abort if users IC violates the periodic BC's. */
    if (init_con[0] != init_con[this->N-1]){
        std::cout << "Error: Chosen initial condition does not satisfy periodic BC!" << std::endl;
        std::cout << "Terminating program..." << std::endl;
    }

    for (int i = 0; i < this->N; i++){
        this->psi_0[i] = init_con[i];
    }
}


void periodic_solver_1d::periodic_euler(){
    // Euler time-stepping in periodic domain
    std::cout << this->dt << std::endl;
}


void periodic_solver_1d::periodic_leapfrog(){
    // Leapfrog time-stepping in periodic domain
}


periodic_solver_1d::~periodic_solver_1d(){
    // No need for cleanup
}
