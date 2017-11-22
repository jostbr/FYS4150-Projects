
#include "basin_solver_1d.hpp"

/* Constructor calling constructor of base class wit input parameters. */
basin_solver_1d::basin_solver_1d(double dx, double dy,
                     int N, double T) : rossby_solver_1d(dx, dy, N, T) {
    // No more initialization in derived class needed
}


/* Function that sets what psi should be at the boundaries. */
void basin_solver_1d::set_boundary_conditions(double bc_0, double bc_N){
    this->bc_0 = bc_0;
    this->bc_N = bc_N;
}


/* Function that sets the initial psi and stores it as a member psi_0. */
void basin_solver_1d::set_initial_condition(double* init_con){
    /* Abort execution if users IC violates users BC's. */
    if (init_con[0] != this->bc_0 || init_con[this->N-1] != this->bc_N){
        std::cout << "Error: Chosen initial condition does not satisfy BC!" << std::endl;
        std::cout << "Terminating program..." << std::endl;
    }

    for (int i = 0; i < this->N; i++){
        this->psi_0[i] = init_con[i];
    }
}


void basin_solver_1d::basin_euler(){
    // Euler time-stepping in basin domain
}


void basin_solver_1d::basin_leapfrog(){
    // Leapfrog time-stepping in basin domain
}


basin_solver_1d::~basin_solver_1d(){
    // No need for cleanup
}
