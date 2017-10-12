#include "solver.h"
#include "planet.h"

#include <vector>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>




solver::solver()
{

}

//Put all classobjects planets into an array
void solver::addPlanet(planet newPlanet)
{
    allPlanets.push_back(newPlanet);
    number_planets += number_planets;

    //return number_planets;
}

void solver::Euler(double Step, double final_time, string filename)
{
    //allPlanets[i].r[0];
    double t = 0.0;
    double r_qubed = 0.0;
    double r_ = 0.0;
    current_r[0] = allPlanets[0].r[0];
    current_r[1] = allPlanets[0].r[1];
    current_r[2] = allPlanets[0].r[2];

    current_v[0] = allPlanets[0].v[0];
    current_v[1] = allPlanets[0].v[1];
    current_v[2] = allPlanets[0].v[2];

    cout << current_r[0] << endl;
    cout << current_v[0] << endl;

    //string filename = output;

    initalize_write_to_file(filename);

    while (t <= final_time){
        next_r[0] = current_r[0] + current_v[0]*Step;
        next_r[1] = current_r[1] + current_v[1]*Step;
        next_r[2] = current_r[2] + current_v[2]*Step;

        r_ = sqrt((current_r[0]*current_r[0])+(current_r[1]*current_r[1])+(current_r[2]*current_r[2]));
        r_qubed = r_*r_*r_;
        cout << "rrr= " << r_qubed<<endl;

        next_v[0] = current_v[0] - (Step*fourpipi*current_r[0])/r_qubed;
        next_v[1] = current_v[1] - (Step*fourpipi*current_r[1])/r_qubed;
        next_v[2] = current_v[2] - (Step*fourpipi*current_r[2])/r_qubed;

        cout << "next_r x= " << next_r[0] <<endl;
        cout << "next_v x= " << next_v[0] <<endl;

        //Need to make a function that writes results to file
        write_row_to_file(t);

        current_r[0] = next_r[0];
        current_r[1] = next_r[1];
        current_r[2] = next_r[2];

        current_v[0] = next_v[0];
        current_v[1] = next_v[1];
        current_v[2] = next_v[2];

        t += Step;


    }

    ofile.close();

}



void solver::write_row_to_file(double t){
    ofile << setw(20) << setprecision(8) << t;
    ofile << setw(20) << setprecision(8) << next_r[0];
    ofile << setw(20) << setprecision(8) << next_r[1];
    ofile << setw(20) << setprecision(8) << next_r[2] << endl;
}


void solver::initalize_write_to_file(string filename){
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "                 t                  x                     y                   z      " << endl;
}

