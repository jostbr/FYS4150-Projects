
#ifndef PERIODIC_SOLVER_1D_HPP
#define PERIODIC_SOLVER_1D_HPP

# include <iostream>
# include <cmath>
# include <armadillo>
# include "rossby_solver.hpp"

class periodic_solver_1d : public rossby_solver {
    public:
        periodic_solver_1d(double dx, int N, double dt, double T, std::string fileout);
        ~periodic_solver_1d();
        void set_initial_condition(double* init_psi, double* init_zeta);
        void periodic_euler();
        void periodic_leapfrog();
};

#endif // PERIODIC_SOLVER_1D_HPP
