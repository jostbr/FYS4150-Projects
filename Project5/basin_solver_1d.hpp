
#ifndef BASIN_SOLVER_1D_HPP
#define BASIN_SOLVER_1D_HPP

# include <iostream>
# include "rossby_solver_1d.hpp"

class basin_solver_1d : public rossby_solver_1d {
    public:
        basin_solver_1d(double dx, double dy, int N, double T);
        ~basin_solver_1d();
        void set_initial_condition(double* init_con);
        void set_boundary_conditions(double bc_0, double bc_N);
        void basin_euler();
        void basin_leapfrog();

    private:
        double bc_0, bc_N;
};

#endif // BASIN_SOLVER_1D_HPP
