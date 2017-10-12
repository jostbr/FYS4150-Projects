//In this main all the functions and constructors for the class is defined
#include "planet.h"

//Defining actions for the empty constructor
planet::planet()
{

    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;

    mass = 1.0;

}
//Constructor with variable values defined in object call?
planet::planet(double x, double y, double z, double vx, double vy, double vz, double M){
    r[0] = x;
    r[1] = y;
    r[2] = z;

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;

    mass = M;
}

//Defining the class functions
double planet::getDistance(planet otherplanet)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->r[0];
    y1 = this->r[1];
    z1 = this->r[2];

    x2 = otherplanet.r[0];
    y2 = otherplanet.r[1];
    z2 = otherplanet.r[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    double rel_distance = sqrt(xx*xx + yy*yy + zz*zz);

    return rel_distance;
}


