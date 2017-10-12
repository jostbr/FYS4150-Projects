
#ifndef PLANET_H
#define PLANET_H

# include <iostream>
# include <cmath>
# include "constants.h"

class planet{
    //static constexpr double mass_sun = 2.0E+30;
    //static constexpr double four_pi_sq = 4*acos(-1.0)*acos(-1.0);

    public:
        planet();
        planet(std::string id, double m, double x, double y, double v_x, double v_y);
        double compute_distance(planet planet_2) const;
        double compute_acc(planet planet_2, int dim) const;
        double compute_total_acc(planet* planets, int num_planets, int dim) const;
        ~planet();

        std::string name;
        double mass;
        double r[2], v[2], a[2];
};

#endif // PLANET_H
