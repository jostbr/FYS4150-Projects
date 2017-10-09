
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>

class planet{
    public:
        planet();
        planet(std::string id, double m, double x, double y, double v_x, double v_y);
        double compute_distance(planet planet_2);
        double compute_force(planet planet_2, int dim);
        ~planet();

        std::string name;
        double mass;
        double mass_sun;
        double r[2], v[2], a[2];
};

#endif // PLANET_H