#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "planet.h"
using namespace std;

//Defining the class solver in this file
class solver
{
    public:

    static constexpr double fourpipi= 4*acos(-1.0)*acos(-1.0);
    static constexpr double G_const = 6.67E-11;     //Nmm/kg

    //To find Center of mass for our system
    double total_mass = 0.0;

    ofstream ofile;

    vector<planet> allPlanets;

    double r_centerofmass[3];
    double x_comp = 0.0;
    double y_comp = 0.0;
    double z_comp = 0.0;

    double sun_velocity[3];



    //Declearing all the classmenmber functions
    solver();

    void addPlanet(planet newPlanet);

    void Euler(double Step, double final_time, string filename);
    void Verlet(double Step, double final_time, string filename);

    void initalize_write_to_file(string filename);
    void write_row_to_file(int i, double t, double x, double y, double z, ofstream **ofiles);

    void resetAcceleration();
    void computeAcceleration();

    void centerofmass();
    void findSolarVelocity();
    void giveSuninitalvelocity();
    void Verlet_CenterofMass(double Step, double final_time, string filename);

};

#endif // SOLVER_H