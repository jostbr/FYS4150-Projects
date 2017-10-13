#ifndef PLANET_H
#define PLANET_H

#include <cmath>
#include <string>

using namespace std;


class planet
{

public:

    static constexpr double fourpipi= 4*acos(-1.0)*acos(-1.0);
    static constexpr double sun_mass = 2.0E30;

    double r[3];
    double v[3];
    double a[3];
    double mass;
    string name;



    //Constructors
    planet();       //Empty constructor
    planet(string id, double x, double y, double z, double vx, double vy, double vz, double M);

    //If I keep class object variables private I need to declear function to access the variables
    //as const //And all functions need to be decleared as friend I think.


    //Class functions - decleartion
    double getDistance(planet otherplanet);
    double Acceleration(planet newPlanet, int dim);
    void resetA();


};




#endif // PLANETS_H
