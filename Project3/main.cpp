
# include <iostream>
# include <cmath>
# include "earth_sun.h"
# include "planet.h"
# include "nbody_solver.h"
# include "constants.h"

int main(){
    //two_body();
    int num_planets = 3;
    planet* planets = new planet[num_planets];

    /* Much used as testing setup. */
    planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 2.0*cnst::pi);
    planet mars("mars", 6.6E+23, 0.0, 1.52, -1.6*cnst::pi, 0.0);
    planet jupiter("jupiter", 1.9E+27, -5.2, 0.0, 0.0, -0.8*cnst::pi);

    /* NASA initial conditions. */
    //planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 2.0*pi);
    //planet mars("mars", 6.6E+23, -1.540775834917753E+00, 6.318374030745555E-01,
    //            -1.7309216762021784, 0.0);
    //planet jupiter("jupiter", 1.9E+27, -5.2, 0.0, 0.0, -0.8*pi);

    planets[0] = earth;
    planets[1] = mars;
    planets[2] = jupiter;

    //std::cout << bodies[1].name << std::endl;

    double t_max = 5.0;       // Upper time in years
    double h = 0.0001;         // Step size in years
    nbody_solver solver(planets, num_planets);
    solver.solve(h, t_max, "verlet");

    delete[] planets;
    return 0;
}
