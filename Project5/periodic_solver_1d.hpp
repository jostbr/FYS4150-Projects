
#ifndef PERIODIC_SOLVER_1D_HPP
#define PERIODIC_SOLVER_1D_HPP

# include <iostream>
# include <cmath>
# include <armadillo>
# include "rossby_solver.hpp"

class periodic_solver_1d : public rossby_solver {
    public:
        periodic_solver_1d(double dx, int N_x, double dt, double T, std::string fileout);
        ~periodic_solver_1d();
        void set_initial_condition(double* init_psi, double* init_zeta);
        void periodic_euler();
        void periodic_leapfrog();

    private:
        void initialize_periodic_matrix(arma::mat& A, int num_rows, int num_cols);
};

#endif // PERIODIC_SOLVER_1D_HPP
