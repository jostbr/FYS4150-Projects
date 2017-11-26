
#ifndef BASIN_SOLVER_1D_HPP
#define BASIN_SOLVER_1D_HPP

# include <iostream>
# include <cmath>
# include "rossby_solver.hpp"
# include "poisson.hpp"

class basin_solver_1d : public rossby_solver {
    public:
        basin_solver_1d(double dx, int N, double dt, double T, std::string fileout);
        ~basin_solver_1d();
        void set_initial_condition(double* init_psi, double* init_zeta);
        void basin_euler();
        void basin_leapfrog();

    private:
        double bc_0, bc_N;
};

#endif // BASIN_SOLVER_1D_HPP
