
#include "rossby_solver_1d.hpp"

rossby_solver_1d::rossby_solver_1d(double dx, double dt, int N, double T){
    this->dx = dx;
    this->dt = dt;
    this->T = T;
    this->N = N;

    alloc_array_1D(this->psi_0, N);
    alloc_array_1D(this->zeta_0, N);
}


/* Function that writes a time and a psi'array to file (on one line). */
void rossby_solver_1d::write_state_to_file(double t, double* psi) const{
    // Write system state to file as row, i.e.: t psi_0, psi_1, ..., psi_N
}

/* Destructor cleaning performing cleanup upon object destruction. */
rossby_solver_1d::~rossby_solver_1d(){
    free_array_1D(this->psi_0);
    free_array_1D(this->zeta_0);
}
