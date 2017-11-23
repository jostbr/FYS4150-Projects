
#ifndef ROSSBY_SOLVERS_1D_HPP
#define ROSSBY_SOLVERS_1D_HPP

# include <iostream>
# include <fstream>
# include <iomanip>
# include "array_alloc.hpp"

class rossby_solver_1d {
    public:
        rossby_solver_1d(double dx, double dt, int N, double T);
        ~rossby_solver_1d();

        double *psi_0, *zeta_0;

    protected:
        void write_state_to_file(double t, double* psi) const;
        double dx, dt, T;
        int N;
};

#endif // ROSSBY_SOLVERS_1D_HPP
