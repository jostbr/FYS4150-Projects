
#include "planet.h"

/* Default constructor function. Initialize everything to zero */
planet::planet(){
    this->name = "";
    this->mass = 0.0;
    this->r[0] = 0.0; this->r[1] = 0.0;
    this->v[0] = 0.0; this->v[1] = 0.0;
    this->a[0] = 0.0; this->a[1] = 0.0;
}

/* Overload constructor function with arguments. */
planet::planet(std::string id, double m, double x, double y, double v_x, double v_y){
    this->name = id;
    this->mass = m;
    this->r[0] = x; this->r[1] = y;
    this->v[0] = v_x; this->v[1] = v_y;
    this->a[0] = 0.0; this->a[1] = 0.0;
}

/* Function that computes distance (sqrt((x_i-x_j)^2 + (y_i-y_j)^2)) between this->planet and planet_2. */
double planet::compute_distance(planet planet_2) const {
    double distance = sqrt((this->r[0] - planet_2.r[0])*(this->r[0] - planet_2.r[0]) +
            (this->r[1] - planet_2.r[1])*(this->r[1] - planet_2.r[1]));
    return distance;
}

/* Function that, through calls to this->compute_distance, computes force acting
 * on this->planet from planet_2 in the dim-direction. */
double planet::compute_acc(planet planet_2, int dim) const {
    double distance_cubed = pow(this->compute_distance(planet_2), 3.0);
    double accel = -((cnst::four_pi_sq*(planet_2.mass/cnst::mass_sun))
                     /distance_cubed)*(this->r[dim] - planet_2.r[dim]);
    return accel;
}

/* Function that, through calls to this->compute_force(), computes total acceleration of
 * this->planet (in dim-direction) due to the sum of gravitational forces from all
 * num_planets planets in array planets*. */
double planet::compute_total_acc(planet* planets, int num_planets, int dim) const {
    planet sun("sun", cnst::mass_sun, 0.0, 0.0, 0.0, 0.0);
    double total_accel = this->compute_acc(sun, dim);     // Include force from sun separately

    for (int i = 0; i < num_planets; i++){
        if (planets[i].name.compare(this->name) != 0){    // No force on itself
            //if (dim == 0) std::cout << "While using force from " << planets[i].name << std::endl;
            total_accel += this->compute_acc(planets[i], dim);
        }
    }

    return total_accel;  // Divide by this->mass to get acceleration
}
/* Note on future efficiency improvement: Compute force from sun separately to save FLOPS. */


/* Destructor function. .*/
planet::~planet(){
    // Destroy stuff
}
