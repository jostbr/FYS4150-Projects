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

    ofstream ofile;

    vector<planet> allPlanets;


    //Declearing all the classmenmber functions
    solver();

    void addPlanet(planet newPlanet);

    void Euler(double Step, double final_time, string filename);
    void Verlet(double Step, double final_time, string filename);

    void initalize_write_to_file(string filename);
    void write_row_to_file(int i, double t, double x, double y, double z, ofstream **ofiles);

    void resetAcceleration();
    void computeAcceleration();
};

#endif // SOLVER_H
