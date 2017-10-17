#include "planet.h"
#include "solver.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <vector>

using namespace std;

int main(int argc, char* argv[]){

    //Read in from terminal arguments
    string filename = argv[1];
    int NumberSteps = atoi(argv[2]);
    float final_time = atof(argv[3]);       //Year

    //Defining Step size
    double Step = final_time/ ((double) NumberSteps);
    cout << "Step= " << Step << endl;

    //I want to crate a class object the Sun - I put it to rest in Origo
    planet Sun("Sun",0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0E30);
    //I want to create a class object planet for Earth- data collected for date 5oct 00.00
    planet Earth("Earth", 9.815271007122528E-01, 2.085477411913488E-01, -1.427693844209846E-04, -3.774829370461290E-03*365, 1.678090463195122E-02*365, -7.297867765979981E-08*365, 6.0E24);


    //To find escape velocity for a two body system
    //planet Earth("Earth", 1.0, 0.0, 0.0, 0.0, 6.3*1.38, 0.0, 6.0E24);

    //Add Jupiter
    planet Jupiter("Jupiter", -4.631702775366028E+00, -2.848650324221476E+00, 1.154126881548682E-01, 3.865209124631586E-03*365, -6.069164573918581E-03*365, -6.123690969326759E-05*365, 1.9E27 );
    //Add Mars
    planet Mars("Mars",-1.505179243871668E+00, 7.131739536964214E-01, 5.169039461400401E-02, -5.425540443407867E-03*365, -1.146705981373801E-02*365, -1.072691469743330E-04*365, 6.39E23 );
    //Add Mercery
    planet Mercury("Mercury", -3.852568085888450E-01, -1.287590592424852E-03, 3.485335827250005E-02, -5.353631439874908E-03*365, -2.691523488282950E-02*365, -1.708985181129974E-03*365, 3.285E23);
    //Add Venus
    planet Venus("Venus", -4.978900478818148E-01, 5.201234856604111E-01, 3.578850004987860E-02, -1.458059295775236E-02, -1.420078774524152E-02, 6.463469308595529E-04, 4.867E24);
    //Add Saturn
    planet Saturn("Saturn", -4.159563835733318E-01, -1.004645250933944E+01, 1.912339898603407E-01, 5.268473041617775E-03, -2.490290405619587E-04, -2.052492138098574E-04, 5.683E26);
    //Add Uranus
    planet Uranus("Uranus", 1.788082074140550E+01, 8.769674540524333E+00, -1.990778755763629E-01, -1.760603538139088E-03, 3.347874640689238E-03, 3.514363762868545E-05,8.681E25);
    //Add Neptun
    planet Neptune("Neptune", 2.860297078564161E+01, -8.857555534925545E+00, -4.767797317610266E-01, 9.080395419979704E-04, 3.017153302252428E-03, -8.344693814272056E-05, 1.024E26);

    solver sys;
    sys.addPlanet(Sun);
//    sys.addPlanet(Mercury);
//    sys.addPlanet(Venus);
    sys.addPlanet(Earth);
//    sys.addPlanet(Mars);
    sys.addPlanet(Jupiter);
//    sys.addPlanet(Saturn);
//    sys.addPlanet(Uranus);
//    sys.addPlanet(Neptune);

    //Solve when we use centerofmass !=Sun
    //Find velocity of SUN that gives total momentum = 0
    sys.findSolarVelocity();
    sys.giveSuninitalvelocity();
    sys.centerofmass();

    sys.Verlet_CenterofMass(Step, final_time, filename);



    //solving Planitary System with Euler
    //sys.Euler(Step, final_time, filename);

    //solving Planitary System with Velocity Verlet
    //sys.Verlet(Step, final_time, filename);



    //Need a destructor off vector allPlanets




};











