
#include "planet.h"

/* Default constructor function. */
planet::planet(){
    name = "";
    mass = 0.0;
    mass_sun = 0.0;
    r[0] = 0.0; r[1] = 0.0;
    v[0] = 0.0; v[1] = 0.0;
    a[0] = 0.0; a[1] = 0.0;
}

/* Overload constructor function. */
planet::planet(std::string id, double m, double x, double y, double v_x, double v_y){
    name = id;
    mass = m;
    mass_sun = 2.0E+30;
    r[0] = x; r[1] = y;
    v[0] = v_x; v[1] = v_y;
    a[0] = 0.0; a[1] = 0.0;
}

double planet::compute_distance(planet planet_2){
    double distance = sqrt((r[0] - planet_2.r[0])*(r[0] - planet_2.r[0]) +
            (r[1] - planet_2.r[0])*(r[1] - planet_2.r[0]));
    return distance;
}

double planet::compute_force(planet planet_2, int dim){
    double four_pi_sq = 4*acos(-1.0)*acos(-1.0);
    double distance = compute_distance(planet_2);
    double force = -four_pi_sq*(planet_2.mass/mass_sun)*mass*(r[dim]
        - planet_2.r[dim])/(pow(distance, 3.0));
    return force;
}

double planet::compute_total_force(planet* planets, int num_planets, int dim){
    double four_pi_sq = 4*acos(-1.0)*acos(-1.0);
    planet sun("sun", 2.0E+30, 0.0, 0.0, 0.0, 0.0);
    double r_sun_cubed = pow(compute_distance(sun), 3.0);
    double total_force = -((four_pi_sq*mass)/(r_sun_cubed))*(r[dim] - sun.r[dim]);

    for (int i = 0; i < num_planets; i++){
        if (planets[i].name.compare(name) != 0){    // No force on itself
            std::cout << "While using force from " << planets[i].name << std::endl;
            total_force += compute_force(planets[i], dim);
        }
    }

    return total_force;
}

/* Destructor function. .*/
planet::~planet(){
    // Destroy stuff
}
