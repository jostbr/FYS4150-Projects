
#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H


class ODE_solver{
    public:
        ODE_solver();
        ~ODE_solver();
    private:
        double h;
        double r[2];
        double v[2];
        double a[2];
};

#endif // ODE_SOLVER_H
