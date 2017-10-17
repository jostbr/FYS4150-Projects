
# include <iostream>
# include <cmath>
# include "earth_sun.h"
# include "planet.h"
# include "nbody_solver.h"
# include "constants.h"

void run_nasa();
void run_three_body();

int main(){
    //two_body();

    run_three_body();

    return 0;
}

void run_nasa(){
    int num_planets = 9;
    planet* planets = new planet[num_planets];

    /* Initial conditions from NASA (JPL) (https://ssd.jpl.nasa.gov/horizons.cgi#results). */
    planet mercury("mercury", 3.3E+23, -3.810263402960284E-01,-1.879843146811533E-01,
                   2.425311737051114, -8.761958449893111);
    planet venus("venus", 4.9E+24, -5.917944681680748E-01, 4.056377351173029E-01,
                 -4.203431900326241, -6.126203994735034);
    planet earth("earth", 6.0E+24, 9.458459982336885E-01, 3.185781851018064E-01,
                 -2.1041651443077454, 5.928179970429281);
    planet mars("mars", 6.6E+23, -1.543025577058266E+00, 6.261351435087505E-01,
                -1.7290323149253264, -4.296563295650305);
    planet jupiter("jupiter", 1.9E+27, -4.606690339852804E+00, 3.185781851018064E-01,
                   1.434541026493422, -2.203707184918641);
    planet saturn("saturn", 5.5E+26, -3.813265876527548E-01, -1.005382158955272E+01,
                  1.9247386951506795, -0.08513862971282533);
    planet uranus("uranus", 8.8E+25, 1.786622934803489E+01, 8.787399621377547E+00,
                  -0.6424714466549416, 1.2191309753735151);
    planet neptune("neptune", 1.03E+26, 2.860706741579854E+01, -8.842134660643014E+00,
                   0.33233181081147617, 1.09939227173294);
    planet pluto("pluto", 1.21E+22, 1.052955154191982E+01, -3.171995053634249E+01,
                 1.1129108016107352, 0.1216767557015243);
    //planet moon("moon", 7.34E+22, 9.453405176570028E-01, 3.209950530345993E-01,
    //            -2.323187760623432, 5.888349836768215);

    planets[0] = mercury;
    planets[1] = venus;
    planets[2] = earth;
    planets[3] = mars;
    planets[4] = jupiter;
    planets[5] = saturn;
    planets[6] = uranus;
    planets[7] = neptune;
    planets[8] = pluto;
    //planets[9] = moon;

    /* User friendly units for parameters. */
    double t_max = 5.0;          // Upper time in years
    double t_write = 10.0;       // Write to file every t_write days
    double h = 1.0/24.0;         // Step size in days

    t_write = t_write/365.0;    // Convert to years before passing argument
    h = h/365.0;                // Convert to years before passing argument

    nbody_solver solar_system(planets, num_planets);      // Create solver object
    solar_system.solve(h, t_max, t_write, "verlet");      // Solve nbody-system using specified method

    delete[] planets;
}

void run_three_body(){
    int num_planets = 4;
    planet* planets = new planet[num_planets];

    /* Much used as testing setup. */
    planet sun("sun", 2.0E+30, 0.0, 0.0, 0.0, 0.0);
    planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 2*cnst::pi);
    planet mars("mars", 6.6E+23, 0.0, 1.52, -1.5*cnst::pi, 0.0);
    planet jupiter("jupiter", 1.9E+27, -5.20, 0.0, 0.0, -0.9*cnst::pi);

    planets[0] = sun;
    planets[1] = earth;
    planets[2] = mars;
    planets[3] = jupiter;

    double t_max = 50.0;          // Upper time in years
    double t_write = 10.0;       // Write to file every t_write days
    double h = 0.0001;         // Step size in days

    t_write = t_write/365.0;    // Convert to years before passing argument
    //h = h/365.0;                // Convert to years before passing argument

    nbody_solver solar_system(planets, num_planets);      // Create solver object
    solar_system.solve(h, t_max, t_write, "verlet");      // Solve nbody-system using specified method

    delete[] planets;
}
