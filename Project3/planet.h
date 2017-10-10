
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
        double compute_total_force(planet* planets, int num_planets, int dim);
        ~planet();

        std::string name;
        double mass;
        double mass_sun;
        double r_curr[2], v_curr[2];
        double r_next[2], v_next[2];
        double a_curr[2];
};

#endif // PLANET_H
