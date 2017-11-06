#include "solver.h"
#include "planet.h"

#include <vector>
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <ctime>



//Empty constructor
solver::solver()
{

}

//Put all classobjects planets into an array allPlanets
void solver::addPlanet(planet newPlanet)
{
    allPlanets.push_back(newPlanet);    //push_back works as append in python
    total_mass += newPlanet.mass;

}

void solver::resetAcceleration(){
    for(int i = 0; i < allPlanets.size(); i++ )
        allPlanets[i].resetA();
}

void solver::computeAcceleration(){
    resetAcceleration();

    for(int i = 0; i < allPlanets.size(); i++){
        for(int j = 0; j < allPlanets.size(); j++){
            if(j == i) continue;
            for(int dim = 0; dim < 3; dim++)
                allPlanets[i].a[dim] += allPlanets[i].Acceleration(allPlanets[j],dim);
        }
    }
}

void solver::Euler(double Step, double final_time, string filename)
{
    double t = 0.0;
    //int print_frame = 0;
    //initalize_write_to_file(filename);
    ofstream* ofiles = new ofstream[allPlanets.size()];
    for (int i = 0; i < allPlanets.size(); i++){
            string filename = allPlanets[i].name;
            filename.append(".txt");
            ofiles[i].open(filename.c_str());
            ofiles[i] << setiosflags(ios::showpoint | ios::uppercase);
            ofiles[i] << "t x y z" << endl;   // Write header to file
            write_row_to_file(i, t, allPlanets[i].r[0], allPlanets[i].r[1], allPlanets[i].r[2], &ofiles);
    }

//    double inital_kin_energy_Earth = 0.5*allPlanets[1].mass*(allPlanets[1].v[0]*allPlanets[1].v[0]
//                            + allPlanets[1].v[1]*allPlanets[1].v[1] + allPlanets[1].v[2]*allPlanets[1].v[2]);

//    double inital_pot_energy_Earth = G_const*((allPlanets[0].mass*allPlanets[1].mass)/allPlanets[0].getDistance(allPlanets[1]));

//    cout << "inital KE = " << inital_kin_energy_Earth << endl;
//    cout << "inital PE = " << inital_pot_energy_Earth << endl;

    while (t <= final_time){
        computeAcceleration();
        t += Step;
        for (int i = 0; i < allPlanets.size(); i++){
            for (int dim=0; dim < 3; dim++){
                allPlanets[i].v[dim] += allPlanets[i].a[dim]*Step;
                allPlanets[i].r[dim] += allPlanets[i].v[dim]*Step;
            }
            write_row_to_file(i,t, allPlanets[i].r[0], allPlanets[i].r[1],allPlanets[i].r[2], &ofiles) ;
        }
    }

//    double final_kin_energy_Earth = 0.5*allPlanets[1].mass*(allPlanets[1].v[0]*allPlanets[1].v[0]
//                               + allPlanets[1].v[1]*allPlanets[1].v[1] + allPlanets[1].v[2]*allPlanets[1].v[2]);
//    double final_pot_energy_Earth = G_const*((allPlanets[0].mass*allPlanets[1].mass)/allPlanets[0].getDistance(allPlanets[1]));

//    cout << "final KE = " << final_kin_energy_Earth << endl;
//    cout << "final PE = " << final_pot_energy_Earth << endl;

//    cout << "Difference Kin energy from inital to final = " << inital_kin_energy_Earth - final_kin_energy_Earth << endl;
//    cout << "Difference Pot energy from inital to final = " << inital_pot_energy_Earth - final_pot_energy_Earth << endl;

//    double inital_tot_energy_Earth = inital_kin_energy_Earth + inital_pot_energy_Earth;
//    double final_tot_energy_Earth = final_kin_energy_Earth + final_pot_energy_Earth;

//    cout << "inital TOTAL energy = " << inital_tot_energy_Earth << endl;
//    cout << "Difference TOTAL energy from inital to final = " << inital_tot_energy_Earth - final_tot_energy_Earth << endl;

    for (int i=0; i<allPlanets.size(); i++){
        ofiles[i].close();
    }
}

//Velocity Verlet
void solver::Verlet(double Step, double final_time, string filename){

    double StepStep = Step*Step;
    double t = 0.0;

    ofstream* ofiles = new ofstream[allPlanets.size()];
    for (int i = 0; i < allPlanets.size(); i++){
            string filename = allPlanets[i].name;
            filename.append(".txt");
            ofiles[i].open(filename.c_str());
            ofiles[i] << setiosflags(ios::showpoint | ios::uppercase);
            ofiles[i] << "t x y z" << endl;   // Write header to file
            write_row_to_file(i, t, allPlanets[i].r[0], allPlanets[i].r[1], allPlanets[i].r[2], &ofiles);
    }

//    //To calculate the inital kinetic energy of the system
//    double total_inital_kin_energy = 0.0;
//    //Starts loop at 1 to avoid including the Sun
//    for (int j=1; j < allPlanets.size(); j++){
//        total_inital_kin_energy +=  0.5*allPlanets[j].mass*(allPlanets[j].v[0]*allPlanets[j].v[0]
//                + allPlanets[j].v[1]*allPlanets[j].v[1] + allPlanets[j].v[2]*allPlanets[j].v[2]);
//    }

//    double total_inital_pot_energy = 0.0;
//    //Can not include the sun as this means deviding on r=0 --> inf
//    for (int j=1; j < allPlanets.size(); j++){
//        total_inital_pot_energy +=  G_const*((allPlanets[0].mass*allPlanets[j].mass)/allPlanets[j].getDistance(allPlanets[0]));
//    }

    while (t <= final_time){
        computeAcceleration();
        t += Step;

        for (int i = 0; i < allPlanets.size(); i++){
            for (int dim=0; dim < 3; dim++){
                if (i==0) continue;         //To keep the Sun fixed in Origo
                //All the current stuff
                allPlanets[i].r[dim] += allPlanets[i].v[dim]*Step + (StepStep/2.0)*allPlanets[i].a[dim];
                //First part of next velocity
                allPlanets[i].v[dim] += allPlanets[i].a[dim]*(Step/2.0);
            }
            computeAcceleration();
            for (int dim=0; dim < 3; dim++){
                if (i==0) continue;         //To keep the Sun fixed in Origo
                allPlanets[i].v[dim] += allPlanets[i].a[dim]*(Step/2.0);
            }
            write_row_to_file(i,t, allPlanets[i].r[0], allPlanets[i].r[1],allPlanets[i].r[2], &ofiles) ;
        }
    }
    for (int i=0; i<allPlanets.size(); i++){
        ofiles[i].close();
    }

//    double total_final_kin_energy = 0.0;
//    //Starts loop at 1 to avoid including the Sun
//    for (int j=1; j < allPlanets.size(); j++){
//        total_final_kin_energy +=  0.5*allPlanets[j].mass*(allPlanets[j].v[0]*allPlanets[j].v[0]
//                + allPlanets[j].v[1]*allPlanets[j].v[1] + allPlanets[j].v[2]*allPlanets[j].v[2]);
//    }


//    cout << "inital KE = " << total_inital_kin_energy << endl;
//    cout << "final KE = " << total_final_kin_energy << endl;

//    //To calculate the potential energy, it should be done with respect to the mass center and total mass?
//    //Under, use sun
//    double total_final_pot_energy = 0.0;
//    //Can not include the sun as this means deviding on r=0 --> inf
//    for (int j=1; j < allPlanets.size(); j++){
//        total_final_pot_energy +=  G_const*((allPlanets[0].mass*allPlanets[j].mass)/allPlanets[j].getDistance(allPlanets[0]));
//    }

//    cout << "inital PE = " << total_inital_pot_energy << endl;
//    cout << "final PE = " << total_final_pot_energy << endl;

//    cout << "inital total energy = " << total_inital_kin_energy+total_inital_pot_energy << endl;
//    cout << "final total energy = " << total_final_kin_energy+total_final_pot_energy << endl;
}


//function to find mass center for system - then put at rest
/* Need to return a vectorial position in r_centerofmass
Includes the sun in the calculations*/
void solver::centerofmass(){
    for (int i = 0; i < allPlanets.size(); i++){
        x_comp += allPlanets[i].r[0]*allPlanets[i].mass;
        y_comp += allPlanets[i].r[1]*allPlanets[i].mass;
        z_comp += allPlanets[i].r[2]*allPlanets[i].mass;
       //r_centerofmass[0]= (allPlanets[i].r[0]*allPlanets[i].mass)/total_mass
    }
    r_centerofmass[0]= x_comp/total_mass;
    r_centerofmass[1]= y_comp/total_mass;
    r_centerofmass[2]= z_comp/total_mass;

    cout << "center off mass in x_dir= " << r_centerofmass[0] << endl;
    cout << "center off mass in y_dir= " << r_centerofmass[1] << endl;
    cout << "center off mass in z_dir= " << r_centerofmass[2] << endl;
}


void solver::findSolarVelocity(){
    double tot_mom_x = 0.0;
    double tot_mom_y = 0.0;
    double tot_mom_z = 0.0;
    //Exlude the sun
    for (int i=1; i< allPlanets.size(); i++){
            tot_mom_x += allPlanets[i].mass*allPlanets[i].v[0];
            tot_mom_y += allPlanets[i].mass*allPlanets[i].v[1];
            tot_mom_z += allPlanets[i].mass*allPlanets[i].v[2];
    }
    //double sun_velocity[3];

    sun_velocity[0] = -tot_mom_x/allPlanets[0].mass;
    sun_velocity[1] = -tot_mom_y/allPlanets[0].mass;
    sun_velocity[2] = -tot_mom_z/allPlanets[0].mass;

    //return sun_velocity;
    cout << "Sun velocity in x dir that make total momentum equal zero = " << sun_velocity[0] << endl;
    cout << "Sun velocity in y dir that make total momentum equal zero = " << sun_velocity[1] << endl;
    cout << "Sun velocity in z dir that make total momentum equal zero = " << sun_velocity[2] << endl;

}

void solver::giveSuninitalvelocity(){
    for (int dim=0; dim<3; dim++){
        allPlanets[0].v[dim]=sun_velocity[dim];
    }
}

//Velocity Verlet -Treates Sun as a regular planet - works for both Sun v_i= 0 and v_i!=0, but sun will move
void solver::Verlet_CenterofMass(double Step, double final_time, string filename){

    double StepStep = Step*Step;
    double t = 0.0;

    ofstream* ofiles = new ofstream[allPlanets.size()];
    for (int i = 0; i < allPlanets.size(); i++){
            string filename = allPlanets[i].name;
            filename.append(".txt");
            ofiles[i].open(filename.c_str());
            ofiles[i] << setiosflags(ios::showpoint | ios::uppercase);
            ofiles[i] << "t x y z" << endl;   // Write header to file
            write_row_to_file(i, t, allPlanets[i].r[0], allPlanets[i].r[1], allPlanets[i].r[2], &ofiles);
    }

//    //To calculate the inital kinetic energy of the system
//    double total_inital_kin_energy = 0.0;
//    //Starts loop at 1 to avoid including the Sun
//    for (int j=1; j < allPlanets.size(); j++){
//        total_inital_kin_energy +=  0.5*allPlanets[j].mass*(allPlanets[j].v[0]*allPlanets[j].v[0]
//                + allPlanets[j].v[1]*allPlanets[j].v[1] + allPlanets[j].v[2]*allPlanets[j].v[2]);
//    }

//    double total_inital_pot_energy = 0.0;
//    //Can not include the sun as this means deviding on r=0 --> inf
//    for (int j=1; j < allPlanets.size(); j++){
//        total_inital_pot_energy +=  G_const*((allPlanets[0].mass*allPlanets[j].mass)/allPlanets[j].getDistance(allPlanets[0]));
//    }

    while (t <= final_time){
        computeAcceleration();
        t += Step;

        for (int i = 0; i < allPlanets.size(); i++){
            for (int dim=0; dim < 3; dim++){
                //if (i==0) continue;         //To keep the Sun fixed in Origo
                //All the current stuff
                allPlanets[i].r[dim] += allPlanets[i].v[dim]*Step + (StepStep/2.0)*allPlanets[i].a[dim];
                //First part of next velocity
                allPlanets[i].v[dim] += allPlanets[i].a[dim]*(Step/2.0);
            }
            computeAcceleration();
            for (int dim=0; dim < 3; dim++){
                //if (i==0) continue;         //To keep the Sun fixed in Origo
                allPlanets[i].v[dim] += allPlanets[i].a[dim]*(Step/2.0);
            }
            write_row_to_file(i,t, allPlanets[i].r[0], allPlanets[i].r[1],allPlanets[i].r[2], &ofiles) ;
        }
    }
    for (int i=0; i<allPlanets.size(); i++){
        ofiles[i].close();
    }

//    double total_final_kin_energy = 0.0;
//    //Starts loop at 1 to avoid including the Sun
//    for (int j=1; j < allPlanets.size(); j++){
//        total_final_kin_energy +=  0.5*allPlanets[j].mass*(allPlanets[j].v[0]*allPlanets[j].v[0]
//                + allPlanets[j].v[1]*allPlanets[j].v[1] + allPlanets[j].v[2]*allPlanets[j].v[2]);
//    }


//    cout << "inital KE = " << total_inital_kin_energy << endl;
//    cout << "final KE = " << total_final_kin_energy << endl;

//    //To calculate the potential energy, it should be done with respect to the mass center and total mass?
//    //Under, use sun
//    double total_final_pot_energy = 0.0;
//    //Can not include the sun as this means deviding on r=0 --> inf
//    for (int j=1; j < allPlanets.size(); j++){
//        total_final_pot_energy +=  G_const*((allPlanets[0].mass*allPlanets[j].mass)/allPlanets[j].getDistance(allPlanets[0]));
//    }

//    cout << "inital PE = " << total_inital_pot_energy << endl;
//    cout << "final PE = " << total_final_pot_energy << endl;
}




void solver::write_row_to_file(int i, double t, double x, double y, double z, ofstream** ofiles){
    (*ofiles)[i] << setw(20) << setprecision(8) << t;
    (*ofiles)[i] << setw(20) << setprecision(8) << x;
    (*ofiles)[i] << setw(20) << setprecision(8) << y;
    (*ofiles)[i] << setw(20) << setprecision(8) << z << endl;
}


void solver::initalize_write_to_file(string filename){
    ofstream* ofiles = new ofstream[allPlanets.size()];
    for (int i = 0; i < allPlanets.size(); i++){
            string filename = allPlanets[i].name;
            filename.append(".txt");
            ofiles[i].open(filename.c_str());
            ofiles[i] << setiosflags(ios::showpoint | ios::uppercase);
            ofiles[i] << "t x y z" << endl;   // Write header to file
            write_row_to_file(i, 0.0, allPlanets[i].r[0], allPlanets[i].r[1], allPlanets[i].r[2], &ofiles);
    }
}









