
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>

class planet{
    static constexpr double mass_sun = 2.0E+30;
    static constexpr double four_pi_sq = 4*acos(-1.0)*acos(-1.0);

    public:
        planet();
        planet(std::string id, double m, double x, double y, double v_x, double v_y);
        double compute_distance(planet planet_2);
        double compute_force(planet planet_2, int dim);
        double compute_acceleration(planet* planets, int num_planets, int dim);
        ~planet();

        std::string name;
        double mass;
        double r[2], v[2];
};

#endif // PLANET_H
