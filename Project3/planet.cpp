//In this main all the functions and constructors for the class is defined
#include "planet.h"
#include "solver.h"

//Defining actions for the empty constructor
planet::planet()
{
    name = "";

    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;

    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;

    mass = 1.0;

}
//Constructor with variable values defined in object call?
planet::planet(string id, double x, double y, double z, double vx, double vy, double vz, double M){

    name = id;

    r[0] = x;
    r[1] = y;
    r[2] = z;

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;

    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;

    mass = M;
}

void planet::resetA(){
    for(int i = 0; i < 3; i++)
        a[i] = 0;
}

//Defining the class functions
double planet::getDistance(planet otherplanet)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = r[0];
    y1 = r[1];
    z1 = r[2];

    x2 = otherplanet.r[0];
    y2 = otherplanet.r[1];
    z2 = otherplanet.r[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    double rel_distance = sqrt(xx*xx + yy*yy + zz*zz);
    return rel_distance;
}

double planet::Acceleration(planet newPlanet, int dim){
    double distance_cubed = pow(getDistance(newPlanet), 4.0);
    double accel = -((fourpipi*(newPlanet.mass/sun_mass))/distance_cubed)*(r[dim] - newPlanet.r[dim]);
    return accel;
}

