
# include <iostream>
# include <cmath>
# include "earth_sun.h"
# include "planet.h"
# include "nbody_solver.h"
# include "constants.h"

/* Function declarations of different setup cases using the model. */
void run_nasa();
void run_three_body();
void run_two_body();

/* Main function calling the various cases. */
int main(){
    //two_body();

    run_two_body();
    //run_three_body();
    //run_nasa();

    return 0;
}

/* Run a simulation of the entire solar system with realistic initial conditions from JPL, NASA. */
void run_nasa(){
    int num_planets = 10;
    planet* planets = new planet[num_planets];

    /* Initial conditions from NASA (JPL) (https://ssd.jpl.nasa.gov/horizons.cgi#results). */
    planet sun("sun", 2.0E+30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    planet mercury("mercury", 3.3E+23, -3.875425742988118E-01, -6.951071548980541E-03, 3.498508347786951E-02,
                   -5.348515291899226E-03*365, -2.692079116485562E-02*365, -1.709103846484009E-03*365);
    planet venus("venus", 4.9E+24, -5.001758135917815E-01, 5.144600047038554E-01, 3.592022525524809E-02,
                 -1.457547680977667E-02*365, -1.420634402726765E-02*365, 6.462282655055163E-04*365);
    planet earth("earth", 6.0E+24, 9.792413350022859E-01, 2.028842602347931E-01, -1.104417905152000E-05,
                 -3.769713222485607E-03*365, 1.677534834992509E-02*365,-1.916440316952949E-07*365);
    planet mars("mars", 6.6E+23, -1.507465009383439E+00, 7.075104721956454E-01, 5.182211933288501E-02,
                -5.420424294970438E-03*365, -1.147261609595497E-02*365, -1.073878125024902E-04*365);
    planet jupiter("jupiter", 1.9E+27, -4.633988541075995E+00, -2.854313805178032E+00, 1.155444133602380E-01,
                   3.870325272607268E-03*365, -6.074720855944709E-03*365, -6.135557504730335E-05*365);
    planet saturn("saturn", 5.5E+26, -4.182421481545518E-01, -1.005211598964372E+01, 1.913657149891211E-01,
                  5.273589189763777E-03*365, -2.545853222131770E-04*365, -2.053678793880784E-04*365);
    planet uranus("uranus", 8.8E+25, 1.787853497569554E+01, 8.764011059567775E+00, -1.989461503709938E-01,
                  -1.755487390163406E-03*365, 3.342318358663112E-03*365, 3.502497227464970E-05*365);
    planet neptune("neptune", 1.03E+26, 2.860068501993164E+01, -8.863219015882104E+00, -4.766480065556578E-01,
                   9.131556899736526E-04*365, 3.011597020226301E-03*365, -8.356560349675644E-05*365);
    planet pluto("pluto", 1.21E+22, 1.050819470206781E+01, -3.172223975950871E+01, 3.537382102084780E-01,
                 3.044230750828683E-03*365, 3.252349625349516E-04*365, -9.049498787998516E-04*365);
    //planet moon("moon", 7.34E+22, 9.817536677155551E-01, 2.029568450722727E-01, -1.516108910221079E-04,
    //            -3.818699428076255E-03, 1.737747002307549E-02, -4.069751782258844E-05);

    planets[0] = mercury;
    planets[1] = venus;
    planets[2] = earth;
    planets[3] = mars;
    planets[4] = jupiter;
    planets[5] = saturn;
    planets[6] = uranus;
    planets[7] = neptune;
    planets[8] = pluto;
    planets[9] = sun;
    //planets[9] = moon;

    /* User friendly units for parameters. */
    double t_max = 10.0;          // Upper time in years
    double t_write = 10.0;       // Write to file every t_write days
    double h = 1.0/24.0;         // Step size in days

    t_write = t_write/365.0;    // Convert to years before passing argument
    h = h/365.0;                // Convert to years before passing argument

    nbody_solver solar_system(planets, num_planets);      // Create solver object
    solar_system.solve(h, t_max, t_write, "verlet");      // Solve nbody-system using specified method

    delete[] planets;
}

void run_three_body(){
    int num_planets = 2;
    planet* planets = new planet[num_planets];

    planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2*cnst::pi, 0.0);
    planet jupiter("jupiter", 1.9E+27, -5.20, 0.0, 0.0, 0.0, -0.9*cnst::pi, 0.0);

    planets[0] = earth;
    planets[1] = jupiter;

    double t_max = 50.0;        // Upper time in years
    double t_write = 4.0;     // Write to file every t_write days
    double h = 1.0/24.0;         // Step size in days

    t_write = t_write/365.0;    // Convert to years before passing argument
    h = h/365.0;                // Convert to years before passing argument

    nbody_solver solar_system(planets, num_planets, true);      // Create solver object
    solar_system.solve(h, t_max, t_write, "verlet");      // Solve nbody-system using specified method

    delete[] planets;
}

void run_two_body(){
    int num_planets = 1;
    planet* planets = new planet[num_planets];

    planet earth("earth", 6.0E+24, 1.0, 0.0, 0.0, 0.0, 2.0*cnst::pi, 0.0);

    planets[0] = earth;

    double t_max = 74.31;        // Upper time in years
    double t_write = 4.0;     // Write to file every t_write days
    double h = 1.0/24.0;         // Step size in days

    t_write = t_write/365.0;    // Convert to years before passing argument
    h = h/365.0;                // Convert to years before passing argument

    nbody_solver solar_system(planets, num_planets, true);      // Create solver object
    solar_system.solve(h, t_max, t_write, "verlet", true);      // Solve nbody-system using specified method

    delete[] planets;
}
