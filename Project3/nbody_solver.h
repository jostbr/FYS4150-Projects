
#ifndef NBODY_SOLVER_H
#define NBODY_SOLVER_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "planet.h"

class nbody_solver{
    public:
        nbody_solver();
        nbody_solver(int body_count, double t_end, double step);
        void euler(planet *nbodies);
        void verlet();
        void write_row_to_file(std::ofstream&, double t, double x, double y);
        ~nbody_solver();

    private:
        int num_bodies;
        double t_max;
        double h;
};

#endif // NBODY_SOLVER_H
