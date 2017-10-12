#include "planet.h"
#include "solver.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>

using namespace std;


/* The Sun and Earth */

int main(int argc, char* argv[]){

    //Read in from terminal arguments
    string filename = argv[1];
    int NumberSteps = atoi(argv[2]);
    float final_time = atof(argv[3]);       //Year

    //Defining Step size
    double Step = final_time/ ((double) NumberSteps);
    cout << "Step= " << Step << endl;

    //I want to create a class object planet for Earth
    planet Earth(9.815271007122528E-01, 2.085477411913488E-01, -1.427693844209846E-04, -3.774829370461290E-03*365, 1.678090463195122E-02*365, -7.297867765979981E-08*365, 6.0E24);
    cout << "inital x position = " << Earth.r[0] << endl; //Do Work

    //I want to crate a class object the Sun - I put it to rest in Origo
    planet Sun(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0E30);

    solver sys;
    sys.addPlanet(Earth);
    sys.addPlanet(Sun);

    //Trying to get the distance r with class function getDistance
    double r = Earth.getDistance(Sun);
    //cout << "r= " << r << endl;   //Do work

    //solving for Earth-Sun with Euler
    sys.Euler(Step, final_time, filename);


}











