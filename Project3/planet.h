
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>
# include "constants.h"

/* Class to represent a celestial body (doesn't have to be a planet). List of members:
 * - name [string] --> One-word description of the planet
 * - mass [double] --> Mass (in kg) of the planet
 * - r[2] [double] --> Position (in AU) in x and y relative to the sun
 * - v[2] [double] --> Velocity (in AU/yr) in x and y
 * - a[2] [double] --> Acceleration (in AU/yr^2) in x and y
 *
 * - compute_distance()  [double] --> Returns distance between wo planets
 * - compute_acc()       [double] --> Returns acceleration to to force from another planet
 * - compute_total_acc() [double] --> Returns net acceleration due to forces from several planets. */

class planet{
    public:
        planet();
        planet(std::string id, double m, double x, double y, double z, double v_x, double v_y, double v_z);
        double compute_distance(planet planet_2) const;
        double compute_acc(planet planet_2, int dim) const;
        double compute_kinetic_energy() const;
        double compute_potential_energy(planet planet_2) const;
        double compute_acc_mercury_GR(planet planet_2, int dim) const;
        ~planet();

        std::string name;
        double mass;
        double r[3], v[3], a[3];
};

#endif // PLANET_H
