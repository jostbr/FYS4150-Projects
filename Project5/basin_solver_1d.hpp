
#ifndef BASIN_SOLVER_1D_HPP
#define BASIN_SOLVER_1D_HPP

# include "rossby_solver_1d.hpp"

class basin_solver_1d : public rossby_solver_1d {
    public:
        basin_solver_1d(double dx, double dy, int N, double T);
        void basin_euler();
        void basin_leapfrog();
};

#endif // BASIN_SOLVER_1D_HPP
