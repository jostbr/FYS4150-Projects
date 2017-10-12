#ifndef PLANET_H
#define PLANET_H

#include <cmath>

using namespace std;


class planet
{

public:
    //double x, y, z, vx, vy, vz, M;

    double r[3];
    double v[3];
    double mass;

    //Constructors
    planet();       //Empty constructor
    planet(double x, double y, double z, double vx, double vy, double vz, double M);

    //If I keep class object variables private I need to declear function to access the variables
    //as const //And all functions need to be decleared as friend I think.


    //Class functions - decleartion
    double getDistance(planet otherplanet);



};




#endif // PLANETS_H
