
#ifndef BASIN_SOLVER_2D_HPP
#define BASIN_SOLVER_2D_HPP

# include <iostream>
# include <cmath>
# include "rossby_solver.hpp"
# include "poisson.hpp"

class basin_solver_2d : public rossby_solver {
    public:
        basin_solver_2d(double dx, double dy, int N_x, int N_y, double dt, double T, std::string fileout);
        ~basin_solver_2d();
        void set_initial_condition(double* init_psi, double* init_zeta);
        void basin_leapfrog();

    private:
        double *bc_0y, *bc_1y, *bc_x0, *bc_x1;
};

#endif // BASIN_SOLVER_2D_HPP
