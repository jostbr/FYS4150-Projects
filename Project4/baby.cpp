#include <iostream>
#include <math.h>       /* exp */
#include <stdlib.h>     /* srand(), rand() */

using namespace std;


//Simple code for solving Ising Model in 2D

int baby()
{
   int L = 2;           //Spins in each dimention
   int N = L*L;         //Total number of Spins

   double start_T = 1.0;    // [k_B T/J]
//   double end_T = 3.0;      // [k_B T/J]
//   double step_T = 0.01;    // Temperature Step

   int MCS = 100;           // Number of MC cycles/sweeps

   // Initalize empty spin matrix
   int spin_matrix[L][L];
   // Fill spin matrix with -1 and 1
   for (int i=0; i<L; i++){
       for (int j=0; j<L; j++){
           int element = rand() % 10 + 0;
           //Can maybe use rand() % -2 +2; to get randum number in interval -1 and 1 -- will not work
           if (element < 5){spin_matrix[i][j] = -1;}
           else {spin_matrix[i][j]=1;}
       cout << spin_matrix[i][j] << endl;
       }
    }

    //Find energy of inital states
    double Energy_init = 0.0;    //Empty varable to fill up with pair energy contributions
    int pbc = L-1;               //Want first matrix element to feel last on same row or column - Periodic Boundary Condition
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            Energy_init += spin_matrix[i][j]*spin_matrix[i][pbc];
            Energy_init += spin_matrix[j][i]*spin_matrix[pbc][i];
            pbc = j;
        }
    }
    cout << "Energy of inital state = " << -Energy_init << endl; //use - because multiply sum by -J

    //Find Magnetization Inital state
    double Mag_init = 0.0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            Mag_init += spin_matrix[i][j];
        }
    }
    cout << "Magnetization of inital state = " << Mag_init << endl;


    //Analytical value of our 2x2 spin system
    //double k_B = 1.38e-23;                  // [J/K]
    double k_B = 1.0;                  // [J/K]
    double beta = 1.0/(k_B*start_T);        //Inverse Temperature
    double Energy_analytical = - (16*exp(8*beta) - 16*exp(-8*beta))/(2*exp(-8*beta) + 2*exp(-8*beta) + 12);
    cout << "Energy by Analytical calculations = " << Energy_analytical << endl;
    //This gives -inf --> probably due to k_B beeing so small that exp --> inf


    //Flip on spin
    int end_index = L;
    int one_index = rand() % end_index + 0;                             //gives random number in interval L-1:0
    cout << "random index = " << one_index << endl;
    int two_index = rand() % end_index + 0;
    cout << "random index = " << two_index << endl;

//    //To test that random indexing actually works
//    for (int i=0; i < 10; i++){
//        cout << "random index = " << rand() % 2 + 0 << endl;            //gives r in interval 1:0 NB
//    }

    cout << "Matrix element before spin flip = " << spin_matrix[one_index][two_index] << endl;
    if (spin_matrix[one_index][two_index] < 0){spin_matrix[one_index][two_index] = 1;}
    else{spin_matrix[one_index][two_index] = -1;}

     cout << "Matrix element after spin flip = " << spin_matrix[one_index][two_index] << endl;

     //Find energy of new states
     double Energy_new = 0.0;    //Empty varable to fill up with pair energy contributions
     //int pbc = L-1;               //Want first matrix element to feel last on same row - Periodic Boundary Condition
     for (int i=0; i<L; i++){
         for (int j=0; j<L; j++){
             Energy_new += spin_matrix[i][j]*spin_matrix[i][pbc];
             Energy_new += spin_matrix[j][i]*spin_matrix[pbc][i];
             pbc = j;
         }

     }
     cout << "Energy of new state = " << -Energy_new << endl;       //uses - because multiply sum with -J
     //But will this work when Energy_new (sum) is negativ also??
     //Seems that way when I test for L=10;



     //Calculate delta_E = Energy_new -  Energy_init
     double delta_E = -Energy_new -  -Energy_init;
     cout << "Energy change by on spin flip = " << delta_E << endl;
     double w;
     int r = rand() % 2 + 0;
     cout << "r = " << r << endl;
     if (delta_E > 0){
         w = exp(-beta*delta_E);
         cout << "w = " << w << endl;
         if (r <= w ){Energy_init = Energy_new;}        //We keep spin configuration

         else{spin_matrix[one_index][two_index] = - spin_matrix[one_index][two_index];} //We change spin configuration back
     }
     //We keep spin configuration
     else{Energy_init = Energy_new;}


     //Find Magnetization Inital state
     double Mag_new = 0.0;
     for (int i=0; i<L; i++){
         for (int j=0; j<L; j++){
             Mag_new += spin_matrix[i][j];
         }
     }
     cout << "Magnetization of new state = " << Mag_new << endl;

     //UPDATE ALL EXPECTATION VALUES AND RUN AGAIN - EACH TIME WE REPEAT IS A MONTE CARLO CYCLE (MCC)
     //AND WE NEED TO DEVIDE OUR EXPECTATION VALUES BY THE NUMBER OF MCC WE PREFORM
     //




    return 0;
}

//Need to include in Qt -
//mpirun -n
//<include> "mpi.h"




//JUST A TEST
//   for (int i=0; i<10; i++){
//       cout << rand() % 10 + 0 << endl;
//   }
