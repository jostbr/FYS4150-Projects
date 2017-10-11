
#ifndef NBODY_SOLVER_H
#define NBODY_SOLVER_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "planet.h"

class nbody_solver{
    public:
        nbody_solver(planet* bodies, int num_bodies);
        void euler(double h, double t_max, std::string method);
        void verlet();
        void initialize_output_files(std::ofstream** ofiles);
        void write_row_to_file(std::ofstream&, double t, double x, double y);
        ~nbody_solver();

    private:
        int num_bodies;
        planet* bodies;
};

#endif // NBODY_SOLVER_H
