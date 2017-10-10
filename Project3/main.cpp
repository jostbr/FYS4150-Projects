
# include <iostream>
# include <cmath>
# include "earth_sun.h"
# include "planet.h"
# include "nbody_solver.h"

int main(){
    //two_body();
    double pi = acos(-1.0);
    static int const num_planets = 1;
    planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 2.0*pi);
    //planet jupiter("jupiter", 1.9E+27, -5.2, 0.0, 0.0, -1.0*pi);
    //planet sun("sun", 2.0E+30, 0.0, 0.0, 0.0, 0.0);

    planet bodies[num_planets];
    bodies[0] = earth;
    //bodies[1] = jupiter;

    //std::cout << bodies[1].name << std::endl;

    double t_max = 2.0;       // Upper time in years
    double h = 0.0001;         // Step size in years
    nbody_solver solver(num_planets, t_max, h);
    solver.euler(bodies);

    return 0;
}
