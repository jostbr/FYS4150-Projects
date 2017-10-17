
#include "planet.h"

/* Default constructor function. Initialize everything to zero */
planet::planet(){
    this->name = "";
    this->mass = 0.0;
    this->r[0] = 0.0; this->r[1] = 0.0; this->r[2] = 0.0;
    this->v[0] = 0.0; this->v[1] = 0.0; this->v[2] = 0.0;
    this->a[0] = 0.0; this->a[1] = 0.0; this->a[2] = 0.0;
}

/* Overload constructor function with arguments. */
planet::planet(std::string id, double m, double x, double y,
               double z, double v_x, double v_y, double v_z){
    this->name = id;
    this->mass = m;
    this->r[0] = x; this->r[1] = y; this->r[2] = z;
    this->v[0] = v_x; this->v[1] = v_y; this->v[2] = v_z;
    this->a[0] = 0.0; this->a[1] = 0.0; this->a[2] = 0.0;
}

/* Function that computes distance (sqrt((x_i-x_j)^2 + (y_i-y_j)^2)) between this->planet and planet_2. */
double planet::compute_distance(planet planet_2) const {
    return sqrt(pow(this->r[0] - planet_2.r[0], 2.0) + pow(this->r[1] - planet_2.r[1], 2.0) +
            pow(this->r[2] - planet_2.r[2], 2.0));
}

/* Function that, through calls to this->compute_distance, computes acceleration for
 * this->planet due to force acting on this->planet from planet_2 in the dim-direction. */
double planet::compute_acc(planet planet_2, int dim) const {
    return -((cnst::four_pi_sq*(planet_2.mass/cnst::mass_sun))
             /pow(this->compute_distance(planet_2), 3.0))*(this->r[dim] - planet_2.r[dim]);
}


/* Destructor function. .*/
planet::~planet(){
    // Destroy stuff
}
