
#ifndef PERIODIC_SOLVER_1D_HPP
#define PERIODIC_SOLVER_1D_HPP

# include <iostream>
# include "rossby_solver_1d.hpp"

class periodic_solver_1d : public rossby_solver_1d {
    public:
        periodic_solver_1d(double dx, double dy, int N, double T, std::string fileout);
        ~periodic_solver_1d();
        void set_initial_condition(double* psi_0);
        void periodic_euler();
        void periodic_leapfrog();
};

#endif // PERIODIC_SOLVER_1D_HPP
