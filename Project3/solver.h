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

    int number_planets = 0;

    //Varaiables needed for Euler
    double current_r[3];
    double current_v[3];
    double next_v[3];
    //Additional Variables needed for velocity Verlet
    double next_r[3];

    ofstream ofile;

    vector<planet> allPlanets;


    //Declearing all the classmenmber functions
    solver();

    void addPlanet(planet newPlanet);

    void Euler(double Step, double final_time, string filename);
    void Verlet(double Step, double final_time, string filename);

    void initalize_write_to_file(string filename);
    void write_row_to_file(double t);

};

#endif // SOLVER_H
