
#ifndef PERIODIC_SOLVER_2D_HPP
#define PERIODIC_SOLVER_2D_HPP

# include <iostream>
# include <cmath>
# include "rossby_solver.hpp"
# include "poisson.hpp"

class periodic_solver_2d : public rossby_solver {
    public:
        periodic_solver_2d(double dx, double dy, int N_x, int N_y, double dt, double T, std::string fileout);
        ~periodic_solver_2d();
        void set_initial_condition(double* init_psi, double* init_zeta);
        void periodic_leapfrog();
};

#endif // PERIODIC_SOLVER_2D_HPP
