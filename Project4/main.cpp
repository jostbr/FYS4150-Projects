#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>       /* exp */
#include <cmath>

#include <ctime>        /* Timing */

#include <random>       /* Random Class */

#include <armadillo>    /* matrix operations */
#include <mpi.h>        /* Parallization */


using namespace std;
using namespace arma;

random_device rd;
mt19937_64 randomEngine(rd());
uniform_real_distribution<double> randomDist(0.0,1.0);



//Declare all the functions
double randomUniform()
{
    return randomDist(randomEngine);
}

void Initialize(int L, mat &Spin_Matrix,  double& Energy, double& Magnetic_Moment);
int PBC(int i, int L, int add);
void WriteResultstoFile(int L, int MCC, double Temperature, vec Expectation_Values);
void Analytical_Values(double Temperature, int L);
void Test_RNG();
void TEST_ALGO(int L, mat &Spin_Matrix_2,  double& Energy_2, double& Magnetic_Moment_2);

// initiate output file
ofstream ofile;

int main(int argc, char* argv[])
{

    int rank, numProcs;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    double start_T;    // [k_B T/J]
    double end_T;      // [k_B T/J]
    double step_T;    // Temperature Step
    int L;

    string filename;
    //Values to be changed Manually between different experiments
    L = 2;
    start_T = 1.0;
    end_T = 1.00001;
    step_T = 0.01;


    if(rank == 0){
        //Test Algorithm
        double Energy_2 = 0.0;
        double Magnetic_Moment_2 = 0.0;
        mat Spin_Matrix_2 = zeros<mat>(2,2);
        cout << "################ Unit Tests ####################################" << endl;
        TEST_ALGO(2, Spin_Matrix_2, Energy_2, Magnetic_Moment_2);

        //Random number generator (RNG) Test
        //Should return a value close to 0.5
        Test_RNG();
    }


    //Initiate empty matrix to hold the expectation values; E, E*E, M, M*M, fabs(M)
    vec Expectation_Values = zeros<mat>(5);             //Locally for each Core
    vec TotExpectation_Values = zeros<mat>(5);          //Store all the Values from all the cores

    if(rank == 0){
        // Declare new file name and add lattice size to file name
        string fileout = filename;
        string argument = to_string(L);
        fileout.append(argument);
        ofile.open(fileout);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "Temperature    Total MCC       E            Cv             M            Chi            |M| " << endl;
    }

   //---------------------------------------------
   //Everything under can go in Metropolis()



    /* Want to time my Algoritm */
    clock_t t;
    if(rank == 0){
        t = clock();
    }



   //Initiate empty variables
   double Energy = 0.0;
   double Magnetic_Moment = 0.0;

   //Initiate zero matrix to fill with spin
   mat Spin_Matrix = zeros<mat>(L,L);

   //Initalize a matrix with random or ground state spin configuration by function
   Initialize(L, Spin_Matrix, Energy, Magnetic_Moment);


   for (double Temperature = start_T; Temperature <= end_T; Temperature += step_T){

       //Values that should be reset for each Temperature
       //int accepted_configurations = 0;
       //int Total_number_MCC = 0;

       //vec Hist_Count = zeros<vec>(1000000);

       //Remove this if-test for large simulations && L=! 2
       if(rank == 0){
           if(L==2){Analytical_Values(Temperature, L);}
       }


       //Empty vector to store the five energy differences
       vec Energy_Change = zeros<vec>(17);
       //Precalculate Energydifferences at a given Temperature
       for (int de=-8; de<=8; de+=4){Energy_Change(de+8)= exp(-de/Temperature);}


       /*To loop over intervall of different Tot_MCC values
         If I only want to use ONE Tot_MCC value;
         Set Tot_MCC = Tot_Tot_MCC */
       int Tot_Tot_MCC = 1250000;

       for (int Tot_MCC = 1250000; Tot_MCC <= Tot_Tot_MCC; Tot_MCC += 250){
           //Loop over all MC cycles

           for (int MCC = 1; MCC <= Tot_MCC; MCC++){
               //Loop over all spins
               for (int i= 0; i<(L*L); i++){
                   //Pick spin at random from Spin_Matrix, floor returns integer value
                   int ix = floor(randomUniform()*L);              //gives random integer number in interval 0:L-1
                   //cout << "random index = " << ix << endl;
                   int iy = floor(randomUniform()*L);
                   //cout << "random index = " << iy << endl;

                   //Find energy of new configuration
                   int Delta_E =  2*Spin_Matrix.at(ix,iy)*
                         (Spin_Matrix.at(ix,PBC(iy,L,-1))
                         +Spin_Matrix.at(PBC(ix,L,-1),iy)
                         +Spin_Matrix.at(ix,PBC(iy,L,1))
                         +Spin_Matrix.at(PBC(ix,L,1),iy));


                   if(randomUniform() <= Energy_Change(Delta_E+8)) {
                       //Flips the spin
                       Spin_Matrix(ix,iy) *= -1;

                       Magnetic_Moment += (double) 2*Spin_Matrix(ix,iy);
                       Energy += (double) Delta_E;
                   }

                   //cout << "Delta_E = " << Delta_E << endl;


               } //Terminates one MC Cycle

           //Want to update all the Expentations Value (five values)
               Expectation_Values(0) += Energy;
               Expectation_Values(1) += Energy*Energy;
               Expectation_Values(2) += Magnetic_Moment;
               Expectation_Values(3) += Magnetic_Moment*Magnetic_Moment;
               Expectation_Values(4) += fabs(Magnetic_Moment);


               //if(rank == 0 && MCC>=5000000){Hist_Count(MCC-5000000)=Energy/L/L;}


           } //Terminates all Tot_MC cycles

//           Total_number_MCC += Tot_MCC;
//           cout << "Total Number of MC sweeps so far = " << Total_number_MCC << endl;

           //Solve part task 4c)
           //cout << "Number of Accepted energy configurations = " << accepted_configurations << endl;

           //cout << "Total Number of MC cycles = " << Tot_MCC << endl;



           //Want to write Expectation values to file
           for(int i = 0; i < 5; i++)
                MPI_Allreduce(&Expectation_Values[i],&TotExpectation_Values[i],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

           if(rank == 0)
                WriteResultstoFile(L, Tot_MCC, Temperature, TotExpectation_Values);

           //----------------------------------------------
           //Metropolis algorithm is ended

           Expectation_Values(0) = 0.0;
           Expectation_Values(1) = 0.0;
           Expectation_Values(2) = 0.0;
           Expectation_Values(3) = 0.0;
           Expectation_Values(4) = 0.0;


       }//Terminates loop over all Tot_Tot_MCC

       /*Uncomment when I want to analyze Probability Distribution
        - Need to comment out the other writing to file code also
       To write results for Hist_Count to file */
//       if (rank==0){
//           string fileout2 = filename;
//           string argument = to_string(L);
//           fileout2.append(argument);
//           ofile.open(fileout2);
//           ofile << setiosflags(ios::showpoint | ios::uppercase);
//           ofile << "Histogram" << endl;
//           for (int i=0; i<(Tot_Tot_MCC); i++){
//               ofile << setw(0) << setprecision(8) << Hist_Count(i) << endl;
//               //cout << Hist_Count(i) << endl;
//           }
//           ofile.close();
//       }



   if(Temperature == end_T && rank == 0){ofile.close();} // close output file}
   }//Terminates for different Temperatures


   if(rank == 0){
       t = clock() - t;
       double time = ((float)t)/CLOCKS_PER_SEC;
       cout << "#################################################" << endl;
       cout << "Time used = " << time << endl;
       cout << "#################################################" << endl;
   }



   MPI_Finalize();




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

  //cout << "Energy of inital state = " << Energy << endl;
  //cout << "Magnetic Moment of inital state = " << Magnetic_Moment << endl;
}

void WriteResultstoFile(int L, int Tot_MCC, double Temperature, vec Expectation_Values)
{
  double norm = 1.0/((double) (4.0*Tot_MCC));  // divided by  number of cycles
  double E = Expectation_Values(0)*norm;
  double EE = Expectation_Values(1)*norm;
  double M = Expectation_Values(2)*norm;
  double MM = Expectation_Values(3)*norm;
  double Mabs = Expectation_Values(4)*norm;
  //cout << setprecision(8) << MM << endl;
  //cout << setprecision(8) << Mabs*Mabs << endl;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double E_variance = (EE- E*E)/L/L;
  double M_variance = (MM- Mabs*Mabs)/L/L;
  //cout << setprecision(8) << M_variance << endl;
  ofile << setw(0) << setprecision(8) << Temperature;
  ofile << setw(15) << setprecision(8) << Tot_MCC*4;
  ofile << setw(15) << setprecision(8) << E/L/L;
  ofile << setw(15) << setprecision(8) << E_variance/Temperature/Temperature;
  ofile << setw(15) << setprecision(8) << M/L/L;
  ofile << setw(15) << setprecision(8) << M_variance/Temperature;
  ofile << setw(15) << setprecision(8) << Mabs/L/L << endl;
} // end output function


void Analytical_Values(double Temperature, int L){

    //The analytical solutions only holds for the 2x2 spin matrix
    double k_B = 1.0;                           // [J/K]
    double beta = 1.0/(k_B*Temperature);        //Inverse Temperature
    double Z_4 = (4*cosh(8*beta)+12);     //*N because we have all the expectation values per spin

    cout << "##########################################################################" << endl;
    cout << "All ANALYTICAL VALUES are per Spin" << endl;
    cout << "##########################################################################" << endl;
    double Energy_analytical =  (-32*sinh(8*beta)/Z_4) / L / L;
    cout << "Energy = " << Energy_analytical << endl;

    double Abs_Mag_Analytical = ((8*exp(8*beta) + 16)/Z_4);
    cout << "Absolut value of magnetic moment  = " << Abs_Mag_Analytical/L/L << endl;

    double Mag_Mag_Analytical = ((32*exp(8*beta)+32)/(Z_4)) ;
    cout << "M*M  = " << Mag_Mag_Analytical/L/L << endl;

    double sigmaE_Analytical = (((256*cosh(8*beta))/Z_4) - ((1024*sinh(8*beta)*sinh(8*beta))/(Z_4*Z_4))) / L /L ;
    cout << "Variance of Energy = " << sigmaE_Analytical<< endl;
    double Cv = ( sigmaE_Analytical  )/(k_B*Temperature*Temperature);
    cout << "Specific Heat  = " << Cv << endl;

    double sigmaM_Analytical = (Mag_Mag_Analytical - Abs_Mag_Analytical*Abs_Mag_Analytical)/ L/L;
    //double sigmaM_Analytical = (((32*exp(8*beta)+32)/Z_4) - ((64*exp(8*beta)*exp(8*beta) + 256)/(Z_4*Z_4))) / L / L ;
    cout << "Variance of Mean magnetization  = " << sigmaM_Analytical << endl;
    double Suseptibility = (  sigmaM_Analytical  )/ (k_B*Temperature);
    cout << "Susceptibility by Analytical calculations = " << Suseptibility << endl;
    cout << "############################################################################" << endl;
}

//I want to make a function that test if RNG is good to go:)
void Test_RNG(){
    int Large_number = 10E6;
    vec Test_numbers = zeros<vec>(Large_number);
    for (int i=0; i<Large_number; i++){
        Test_numbers(i) = randomUniform();
    }


    //Want to check if mean is 1/2 for uniform distribution
    double Test_result = mean(Test_numbers);
    cout << "Mean of random numbers = " << Test_result << endl;
    double allowed_offset = 0.01;
    if(fabs(Test_result - 0.5) < allowed_offset){
        cout << "Random numer generator: OK" << endl;
    }
    else{
        cout << "Random numer generator: NOT OK" << endl;
        exit(EXIT_FAILURE);
    }


    //To proper test RNG, increase Large_number and repeat many times, then mean of all the means:)
}

void TEST_ALGO(int L, mat &Spin_Matrix_2,  double& Energy_2, double& Magnetic_Moment_2)
{
    // Setup spin matrix for GROUND STATE T < 1.5
    for(int x =0; x < L; x++) {
      for (int y= 0; y < L; y++){
          Spin_Matrix_2(x,y) = 1;
          //Setup initial magnetization
          Magnetic_Moment_2 +=  (double) Spin_Matrix_2(x,y);
      }
    }
    // Setup initial energy
    for(int x = 0; x < L; x++) {
      for (int y = 0; y < L; y++){
        Energy_2 -= (double) Spin_Matrix_2(x,y)*(Spin_Matrix_2(PBC(x,L,-1),y) + Spin_Matrix_2(x,PBC(y,L,-1)));
      }
    }

    //Ground state - All spins points up
    //TEST_ALGO should return Energy_2 = -8
    //TEST_ALGO should return Magnetic_Moment_2 = 4
    if(Energy_2 == -8.0 & Magnetic_Moment_2 == + 4.0){
        cout << "Test passed for a 2x2 spin matrix in Ground State" << endl;
    }
    else{
        cout << "Test NOT passed for a 2x2 spin matrix in Ground State" << endl;
        exit(EXIT_FAILURE);
    }

}


