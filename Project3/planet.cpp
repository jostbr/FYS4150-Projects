
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
             /pow(this->compute_distance(planet_2), 3.5))*(this->r[dim] - planet_2.r[dim]);
}

/* Function that, through calls to this->compute_distance, computes acceleration for
 * mercury taking general relativity effects into account from planet_2 in dim-direction. */
double planet::compute_acc_mercury_GR(planet planet_2, int dim) const {
    double l_x = (this->r[1]*this->v[2] - this->r[2]*this->v[1]);
    double l_y = (this->r[2]*this->v[0] - this->r[0]*this->v[2]);
    double l_z = (this->r[0]*this->v[1] - this->r[1]*this->v[0]);
    double l = sqrt(l_x*l_x + l_y*l_y + l_z*l_z);
    double r = this->compute_distance(planet_2);

    double acc = -(((cnst::four_pi_sq*(planet_2.mass/cnst::mass_sun))/pow(r, 3.0))*
                   (this->r[dim] - planet_2.r[dim]))*(1.0 + 3.0*l*l/(r*r*cnst::c*cnst::c));
    return acc;
}

/* Function that returns the kinetic energy of the planet K_E = 0.5*m*v^2. */
double planet::compute_kinetic_energy() const {
    return 0.5*this->mass*(pow(this->v[0], 2.0) + pow(this->v[1], 2.0) + pow(this->v[2], 2.0));
}

double planet::compute_potential_energy(planet planet_2) const {
    return -cnst::G*this->mass*planet_2.mass/this->compute_distance(planet_2);
}


/* Destructor function. .*/
planet::~planet(){
    // Destroy stuff
}
