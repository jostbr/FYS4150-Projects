#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>       /* exp */
#include <cmath>


#include <random>

#include <armadillo>    /* matrix operations */



using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> randomDist(0.0,1.0);



//Declare all the functions
double randomUniform()
{
    return randomDist(randomEngine);
}

void Initialize(int L, mat &Spin_Matrix,  double& Energy, double& Magnetic_Moment);
int PBC(int i, int L, int add);
void WriteResultstoFile(int L, int MCC, double Temperature, vec Expectation_Values);
void Analytical_Values(double Temperature, int L, int Tot_MCC);
void Test_RNG();

// initiate output file
ofstream ofile;

int main(int argc, char* argv[])
{
    double start_T;    // [k_B T/J]
    double end_T;      // [k_B T/J]
    double step_T;    // Temperature Step
    int L;

    int Tot_MCC;           // Number of MC cycles/sweeps

    string filename=argv[1];
    L = atoi(argv[2]);
    Tot_MCC = atoi(argv[3]);
    start_T = atof(argv[4]);
    end_T = atof(argv[5]);
    step_T = atof(argv[6]);

    //Test_RNG();


   //int N = L*L;         //Total number of Spins


    //Initiate empty matrix to hold the expectation values; E, E*E, M, M*M, fabs(M)
    vec Expectation_Values = zeros<mat>(5);

    // Declare new file name and add lattice size to file name
    string fileout = filename;
    string argument = to_string(L);
    fileout.append(argument);
    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "Temperature    Tot_MCC       E            Cv             M            Chi            |M| " << endl;


    //---------------------------------------------
   //Everything under can go in Metropolis() - work in progress



   //Initiate empty vector to hold the  for differnt energy changes (de)
   //vec Energy_Bin = zeros<vec>(L*L*4);

   //Initiate empty variables
   double Energy = 0.0;
   double Magnetic_Moment = 0.0;
   //Initiate zero matrix to fill with spin
   mat Spin_Matrix = zeros<mat>(L,L);

   //Initalize a matrix with random or ground state spin configuration by function
   Initialize(L, Spin_Matrix, Energy, Magnetic_Moment);

   //Precalculate Energydifferences at a given Temperature
   for (double Temperature = start_T; Temperature <= end_T; Temperature += step_T){

       //Values that should be reset for each Temperature
       int accepted_configurations = 0;
       int Total_number_MCC = 0;

       //Remove this if-test for large simulations..
       if(L==2){Analytical_Values(Temperature, L, Tot_MCC);}

       //Empty vector to store the five energy differences
       vec Energy_Change = zeros<vec>(17);
       for (int de=-8; de<=8; de+=4){Energy_Change(de+8)= exp(-de/Temperature);}
       //cout << Energy_Change << endl;



       int Tot_Tot_MCC = 10000000;
       for (int Tot_MCC = 10000000; Tot_MCC <= Tot_Tot_MCC; Tot_MCC += 100000){
           //Loop over all MC cycles
           for (int MCC = 0; MCC <= Tot_MCC; MCC++){
               //Loop over all spins
               for (int i= 0; i<(L*L); i++){
                   //Pick spin at random from Spin_Matrix
                   int ix = floor(randomUniform()*L);                             //gives random number in interval L-1:0
                   //cout << "random index = " << ix << endl;
                   int iy = floor(randomUniform()*L);
                   //cout << "random index = " << iy << endl;

                   Spin_Matrix(ix,iy) *= -1;    //Flips the spin we chose at random

                   //Find energy of new configuration
                   int Delta_E =  -2*Spin_Matrix.at(ix,iy)*
                         (Spin_Matrix.at(ix,PBC(iy,L,-1))+Spin_Matrix.at(PBC(ix,L,-1),iy)+Spin_Matrix.at(ix,PBC(iy,L,1))+Spin_Matrix.at(PBC(ix,L,1),iy));

                   //cout << "Delta_E = " << Delta_E << endl;

                   //Want to change new spin configuration back if condition under is met
                   if (Energy_Change(Delta_E+8) <= randomUniform() ){
                       Spin_Matrix(ix,iy) *= -1;

                   }
                   else{
                       Energy += (double) Delta_E;
                       Magnetic_Moment += (double) 2*Spin_Matrix(ix,iy);

                       accepted_configurations +=1;

//                       if(Delta_E == -8){ Bin(0) +=1;}
//                       if(Delta_E == -4){ Bin(1) +=1;}
//                       if(Delta_E == 0){ Bin(2) +=1;}
//                       if(Delta_E == 4){ Bin(3) +=1;}
//                       if(Delta_E == 8){ Bin(4) +=1;}
                   }


               } //Terminates one MC Cycle

           //Want to update all the Expentations Value (five values)
               Expectation_Values(0) += Energy;
               Expectation_Values(1) += Energy*Energy;
               Expectation_Values(2) += Magnetic_Moment;
               Expectation_Values(3) += Magnetic_Moment*Magnetic_Moment;
               Expectation_Values(4) += fabs(Magnetic_Moment);

               //Energy_Bin(Energy+(L*L*2)) += 1;



           } //Terminates all Tot_MC cycles

           Total_number_MCC += Tot_MCC;
           cout << "Total Number of MC sweeps so far = " << Total_number_MCC << endl;

           //Solve part task 4c)
           cout << "Number of Accepted energy configurations = " << accepted_configurations << endl;

           //Solve task 4d) HISTOGRAM as function of T=1.0 and T=2.4 (and tot number of MCC)
//           cout << "Number of bin elements for DeltaE -8 = " << Bin(0) << endl;
//           cout << "Number of bin elements for DeltaE -4 = " << Bin(1) << endl;
//           cout << "Number of bin elements for DeltaE 0 = " << Bin(2) << endl;
//           cout << "Number of bin elements for DeltaE 4 = " << Bin(3) << endl;
//           cout << "Number of bin elements for DeltaE 8 = " << Bin(4) << endl;

           //cout << "Total Number of MC cycles = " << Tot_MCC << endl;

           //----------------------------------------------

           //Want to write Expectation values to file
           WriteResultstoFile(L, Tot_MCC, Temperature, Expectation_Values);

           Expectation_Values(0) = 0.0;
           Expectation_Values(1) = 0.0;
           Expectation_Values(2) = 0.0;
           Expectation_Values(3) = 0.0;
           Expectation_Values(4) = 0.0;

       }//Terminates loop over all Tot_Tot_MCC

//       for (int i=0; i<(L*L); i++){
//           cout << Energy_Bin(i) << endl;
//       }

//       Expectation_Values(0) = 0.0;
//       Expectation_Values(1) = 0.0;
//       Expectation_Values(2) = 0.0;
//       Expectation_Values(3) = 0.0;
//       Expectation_Values(4) = 0.0;

   if(Temperature == end_T){ofile.close();} // close output file}
   }//Terminates for different Temperatures


   return 0;
}




// Function for PeriodicBoundary conditions
int PBC(int i, int L, int add)
{
  return (i+L+add) % (L);   //This returns the reminder of the division (i+L+add)/L
}

//Function to initialize a spin configuration in a Spin_Matrix
void Initialize(int L, mat &Spin_Matrix,  double& Energy, double& Magnetic_Moment)
{
  // Setup random spin matrix
//  for(int x =0; x < L; x++) {
//    for (int y= 0; y < L; y++){
//        double element = randomUniform();
//        //cout << element << endl;
//        if (element < 0.5){Spin_Matrix(x,y) = -1;}
//        else {Spin_Matrix(x,y) = 1;}

//        //cout << Spin_Matrix(x,y) << endl;

//        //Setup initial magnetization
//        Magnetic_Moment +=  (double) Spin_Matrix(x,y);
//    }
//  }

  // Setup spin matrix for GROUND STATE T < 1.5
  for(int x =0; x < L; x++) {
    for (int y= 0; y < L; y++){
        Spin_Matrix(x,y) = 1;

        //Setup initial magnetization
        Magnetic_Moment +=  (double) Spin_Matrix(x,y);
    }
  }


  // Setup initial energy
  for(int x = 0; x < L; x++) {
    for (int y = 0; y < L; y++){
      Energy -= (double) Spin_Matrix(x,y)*(Spin_Matrix(PBC(x,L,-1),y) + Spin_Matrix(x,PBC(y,L,-1)));
    }
  }

  cout << "Energy of inital state = " << Energy << endl;
  cout << "Magnetic Moment of inital state = " << Magnetic_Moment << endl;
}

void WriteResultstoFile(int L, int Tot_MCC, double Temperature, vec Expectation_Values)
{
  double norm = 1.0/((double) (Tot_MCC));  // divided by  number of cycles
  double E = Expectation_Values(0)*norm;
  double EE = Expectation_Values(1)*norm;
  double M = Expectation_Values(2)*norm;
  double MM = Expectation_Values(3)*norm;
  double Mabs = Expectation_Values(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double E_variance = (EE- E*E)/L/L;
  double M_variance = (MM - Mabs*Mabs)/L/L;
  ofile << setw(0) << setprecision(8) << Temperature;
  ofile << setw(15) << setprecision(8) << Tot_MCC;
  ofile << setw(15) << setprecision(8) << E/L/L;
  ofile << setw(15) << setprecision(8) << E_variance/Temperature/Temperature;
  ofile << setw(15) << setprecision(8) << M/L/L;
  ofile << setw(15) << setprecision(8) << M_variance/Temperature;
  ofile << setw(15) << setprecision(8) << Mabs/L/L << endl;
} // end output function


void Analytical_Values(double Temperature, int L, int Tot_MCC){
    //double norm = 1.0/((double) (Tot_MCC));  // divided by  number of cycles

    //The analytical solutions only holds for the 2x2 spin matrix
    double k_B = 1.0;                           // [J/K]
    double beta = 1.0/(k_B*Temperature);        //Inverse Temperature
    double Z_4 = (4*cosh(8*beta)+12);     //*N because we have all the expectation values per spin

    double Energy_analytical =  (-32*sinh(8*beta)/Z_4) / L / L;
    cout << "Energy by Analytical calculations = " << Energy_analytical << endl;

    double Abs_Mag_Analytical = ((8*exp(8*beta) + 16)/Z_4);
    cout << "Absolut value of magnetic moment by Analytical calculations = " << Abs_Mag_Analytical/L/L << endl;

    double Mag_Mag_Analytical = ((32*exp(8*beta)+32)/(Z_4)) ;
    cout << "M*M by Analytical Calculation = " << Mag_Mag_Analytical/L/L << endl;

    double sigmaE_Analytical = (((256*cosh(8*beta))/Z_4) - ((1024*sinh(8*beta)*sinh(8*beta))/(Z_4*Z_4))) / L /L ;
    cout << "Variance of Energy by Analytical calculations = " << sigmaE_Analytical<< endl;
    double Cv = ( sigmaE_Analytical  )/(k_B*Temperature*Temperature);
    cout << "Specific Heat by Analytical calculations = " << Cv << endl;

    double sigmaM_Analytical = (Mag_Mag_Analytical - Abs_Mag_Analytical*Abs_Mag_Analytical)/ L/L;
    //double sigmaM_Analytical = (((32*exp(8*beta)+32)/Z_4) - ((64*exp(8*beta)*exp(8*beta) + 256)/(Z_4*Z_4))) / L / L ;
    cout << "Variance of Mean magnetization by Analytical calculations = " << sigmaM_Analytical << endl;
    double Suseptibility = (  sigmaM_Analytical  )/ (k_B*Temperature);
    cout << "Susceptibility by Analytical calculations = " << Suseptibility << endl;
}

//I want to make a function that test if RNG is good to go:)

void Test_RNG(){
    int Large_number = 10E1;
    vec Test_numbers = zeros<vec>(Large_number);
    for (int i=0; i<=Large_number; i++){
        Test_numbers(i) = randomUniform();
    }


    //Want to check if integral gives 1/2
    //vec X = linspace<vec>(0, 1.0, 1600);

   // mat Z =  trapz(X,Test_numbers);
    //double Z = trapz(Test_numbers);

   // cout << "Integral of random numbers = " << Z << endl;

}



