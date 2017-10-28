
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>
# include "constants.h"

/* Class to represent a celestial body (doesn't have to be a planet). List of members:
 * - name [string] --> One-word description of the planet
 * - mass [double] --> Mass (in kg) of the planet
 * - r[3] [double] --> Position (in AU) in x, y and z relative to origin
 * - v[3] [double] --> Velocity (in AU/yr) in x, y and z
 * - a[3] [double] --> Acceleration (in AU/yr^2) in x, y and z
 *
 * - compute_distance()         [double] --> Returns distance between wo planets
 * - compute_acc()              [double] --> Returns acceleration to to force from another planet
 * - compute_kinetic_energy()   [double] --> Returns kinetic energy of the planet
 * - compute_potential_energy() [double] --> Returns potential energy of a planet relative to another planet
 * - compute_acc_mercury_GR()   [double] --> Returns acceleration of planet taking GR effects into account */

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
