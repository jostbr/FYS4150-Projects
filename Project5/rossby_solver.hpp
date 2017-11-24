
#ifndef ROSSBY_SOLVERS_1D_HPP
#define ROSSBY_SOLVERS_1D_HPP

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "array_alloc.hpp"

class rossby_solver {
    public:
        rossby_solver(double dx, int N, double dt, double T, std::string fileout);
        ~rossby_solver();
        void get_initial_condition(double* psi_holder, double* zeta_holder) const;
        void display_model_config() const;

    protected:
        void write_state_to_file(double t, double* psi);

        double *psi_0, *zeta_0;
        double dx, dt, T;
        int N;

        std::string fileout;
        std::ofstream outfile;
};

#endif // ROSSBY_SOLVERS_1D_HPP
