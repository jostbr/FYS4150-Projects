
#include "periodic_solver_1d.hpp"

periodic_solver_1d::periodic_solver_1d(double dx, double dy,
                        int N, double T) : rossby_solver_1d(dx, dy, N, T) {
    // No more initialization in derived class needed
}

void periodic_solver_1d::periodic_euler(){
    // Euler time-stepping in periodic domain
}

void periodic_solver_1d::periodic_leapfrog(){
    // Leapfrog time-stepping in periodic domain
}
