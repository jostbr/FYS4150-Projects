
#ifndef NBODY_SOLVER_H
#define NBODY_SOLVER_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <ctime>
# include "planet.h"
# include "constants.h"

class nbody_solver{
    public:
        nbody_solver(planet* bodies, int num_bodies);
        void solve(double h, double t_max, double t_write, std::string method);
        void write_row_to_file(int file_index, double t, double x, double y);
        ~nbody_solver();

    private:
        void euler(double h, double t_max, int frame_write);
        void verlet(double h, double t_max, int frame_write);
        double compute_total_acc(planet subject, planet* objects, int dim) const;
        int num_bodies;
        planet* bodies;
        std::ofstream* ofiles;
};

#endif // NBODY_SOLVER_H
