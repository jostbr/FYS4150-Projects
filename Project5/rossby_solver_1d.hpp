
#ifndef ROSSBY_SOLVERS_1D_HPP
#define ROSSBY_SOLVERS_1D_HPP


class rossby_solver_1d {
    public:
        rossby_solver_1d(double dx, double dt, int N, double T);
        ~rossby_solver_1d();

    protected:
        void write_state_to_file(double t, double* psi) const;
        double dx, dt, T;
        int N;

};

#endif // ROSSBY_SOLVERS_1D_HPP
