
#include "basin_solver_1d.hpp"

basin_solver_1d::basin_solver_1d(double dx, double dy,
                     int N, double T) : rossby_solver_1d(dx, dy, N, T) {
    // No more initialization in derived class needed
}

void basin_solver_1d::basin_euler(){
    // Euler time-stepping in basin domain
}

void basin_solver_1d::basin_leapfrog(){
    // Leapfrog time-stepping in basin domain
}
