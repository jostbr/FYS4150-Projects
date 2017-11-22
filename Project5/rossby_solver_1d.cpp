
#include "rossby_solver_1d.hpp"

rossby_solver_1d::rossby_solver_1d(double dx, double dt, int N, double T){
    this->dx = dx;
    this->dt = dt;
    this->T = T;
    this->N = N;
}

void rossby_solver_1d::write_state_to_file(double t, double* psi) const{
    // Write system state to file as row, i.e.: t psi_0, psi_1, ..., psi_N
}

rossby_solver_1d::~rossby_solver_1d(){
    // No cleanup needed
}
