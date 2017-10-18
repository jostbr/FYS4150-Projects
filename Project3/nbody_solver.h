
#ifndef NBODY_SOLVER_H
#define NBODY_SOLVER_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <cmath>
# include <ctime>
# include "planet.h"
# include "constants.h"

class nbody_solver{
    public:
        nbody_solver(planet* bodies, int n, bool implicit_sun = false);
        void solve(double h, double t_max, double t_write, std::string method);
        void write_row_to_file(int file_index, double t, double x, double y, double z);
        ~nbody_solver();

    private:
        void euler(double h, double t_max, int frame_write);
        void verlet(double h, double t_max, int frame_write);
        double compute_total_acc(planet subject, planet* objects, int dim) const;
        void display_kinetic_energy(double time) const;
        int num_bodies;
        bool fixed_sun;
        planet* bodies;
        std::ofstream* ofiles;
};

#endif // NBODY_SOLVER_H
