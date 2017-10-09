
#include "planet.h"

/* Default constructor function. */
Planet::Planet(){
    planet_name = "";
    r[0] = 0.0; r[1] = 0.0;
    v[0] = 0.0; v[1] = 0.0;
    a[0] = 0.0; a[1] = 0.0;
}

/* Overload constructor function. */
Planet::Planet(std::string planet_name, double x, double y, double v_x, double v_y){
    planet_name = "";
    r[0] = x; r[1] = y;
    v[0] = v_x; v[1] = v_y;
    a[0] = 0.0; a[1] = 0.0;
}

/* Destructor function. .*/
Planet::~Planet(){
    // Destroy stuff
}
