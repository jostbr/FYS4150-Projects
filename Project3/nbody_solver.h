
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

/* Class to represent an n-body system and its solutions. List of members:
 * - num_bodies [int] --> Number of bodies in simulation
 * - fixed_sun  [bool] --> Whether to use fixed sun or moveable sun
 * - bodies     [planet*] --> Array of planets
 * - ofiles     [doublesrd::ofstream*] --> Array of file objects for output files
 *
 * - solve()                  [void] --> MAnager function for solving n-body system. Also performs tests
 * - write_row_to_file()      [void] --> Writes a row of t, x, y, z values to output file
 * - euler                    [void] --> Euler algorithm for time stepping the n-body system
 * - verlet()                 [void] --> Velocity Verlet algorithm for time stepping the n-body system
 * - compute_total_acc        [double] --> Returns total acceleration on one body
 * - compute_energy           [void] --> Prints out the total kinetic, potential and mechanical energy
 * - compute_angular_momentum [void] --> Prints out the total angular momentum of the system
 * - compute_center_mass      [void] --> Prints out coordinates of the center of mass for the system */

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
        void compute_energy(double time) const;
        void compute_angular_momentum(double time) const;
        void compute_center_mass() const;
        int num_bodies;
        bool fixed_sun;
        planet* bodies;
        std::ofstream* ofiles;
};

#endif // NBODY_SOLVER_H
