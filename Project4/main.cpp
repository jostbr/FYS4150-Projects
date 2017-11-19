
# include <string>
# include <mpi.h>
# include "ising.h"
# include "mem_alloc.h"
# include "unit_tests.h"

int main(int argc, char* argv[]){
    int my_rank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_rank == 0){
        std::cout << "\nRUNNING TESTS..." << std::endl;
        std::cout << "===============================================" << std::endl;
        TEST_get_periodic_index();
        TEST_compute_energy_and_moment();
        TEST_initialize_spin_config_ordered();
        TEST_initialize_spin_config_rng();
        std::cout << "===============================================" << std::endl;

        std::cout << "\nAnalytical values for 2x2 lattice (per spin) with T = 1.0" << std::endl;
        std::cout << "----------------------------------------------------------------" << std::endl;
        print_analytical_values(1.0);
    }

    double T_0 = 2.2;
    double T_N = 2.35;
    double dT = 0.02;
    int num_spins = 80;
    int num_mc_cycles = 4000000;

    if (num_mc_cycles % num_procs != 0){    // If num_mc_cycles can't be evenly distributed
        std::cout << "Error: Choose num_mc_cyles and num_procs --> num_mc_cycles % num_procs = 0" << std::endl;
        std::cout << "Terminating program...\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string fileout_01 = "PR_80x80_TR.txt";

    ising_mpi(T_0, T_N, dT, num_spins, num_mc_cycles, my_rank, num_procs, fileout_01);

    MPI_Finalize();
}

/* Function implementing a parallel (MPI) Ising model performing Monte Carlo simulations. Takes input
 * initial and final temperature of the system, a temperature step, number of spins in one direction (of
 * the 2D square lattice) and the number of Monte Carlo cycle one wishes to perform (the more the better).
 * Also, for the parallelization, it needs the rank of the caller and the total number of processes. */
void ising_mpi(double T_min, double T_max, double dT, int num_spins, int num_mc_cycles,
               int my_rank, int num_procs, std::string fileout){
    /* ================================== Variable declarations ================================== */
    double E, M;                    // To hold energy and magentic moment each temperature iteration
    double T = T_min;               // Temperature starts at minimum
    long long accepted_states = 0;        // To count number of total accepted states
    double averages[5];             // Array to hold total expectation values of E, E^2, M, M^2, |M|
    double exp_dE[5];               // Array to hold pre-calculated w = exp(-dE/T)
    int delta_E[5];                 // Array to hold pre-caclualted energy differences
    int** spins;                    // Spin array (2D) to hold +1 or -1 values
    int total_num_spins = num_spins*num_spins;

    /* Variables especially important on a per process-level. */
    double my_averages[5];              // Array to hold local expectation values of E, E^2, M, M^2, |M|
    long long my_accepted_states;             // To count local number of total accepted states
    int my_num_mc_cycles = num_mc_cycles/num_procs;     // Assume integer division gives zero remainder
    std::cout << "Number of MC cycles performed by process "
              << my_rank << ": " << my_num_mc_cycles << std::endl;

    /* ================================== Pre-loop setup ================================== */
    alloc_array_2D(spins, num_spins, num_spins);    // Allocate memory for 2D spin array.

    delta_E[0] = -8; delta_E[1] = -4; delta_E[2] = -0; delta_E[3] = 4; delta_E[4] = 8;  // Possible dE's

    std::ofstream outfile, outfile_02;  // Only used by process 0

    /* Open output file and write header info to file. */
    if (my_rank == 0){
        outfile.open(fileout);
        outfile_02.open("energy_count.txt");
        outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        outfile_02 << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        outfile << std::setw(10) << "Num spins" << std::setw(18) << "Tempreature" << std::setw(18)
                << "MC cycles" << std::setw(18) << "Mean energy" << std::setw(18) << "Mean abs mag"
                << std::setw(18) << "Mean mag" << std::setw(18) << "Energy var" << std::setw(18)
                << "Mean abs var" << std::setw(18) << "Heat cap" << std::setw(18)
                << "Mag susc" << std::setw(18) << "Acc states"<< std::endl;
    }

    std::random_device rd;      // Instantiate random device
    std::mt19937_64 rng(rd());  // Instantiate random number generator
    std::uniform_real_distribution<double> uniform(0.0, 1.0);  // Use a uniform dist for the generator

    initialize_spin_config_rng(spins, num_spins);       // Set initial spin config for current T

    double time_end, time_used, time_start = MPI_Wtime();   // Start timing on all processes

    /* Main loop over temperature doing Monte Carlo simulations each temperature step. */
    /* ================================================================================= */
    while (T <= T_max){
        std::cout << T << std::endl;
        /* ========================== Pre-MC-cycles initialization =========================== */
        E = 0;      // Start energy at zero for every temperature
        M = 0;      // Start magnetic moment at zero for every temperature
        my_accepted_states = 0;   // Reset number of accepted states every temperature

        for (int m = 0; m < 5; m++){
            exp_dE[m] = exp(-delta_E[m]/T);     // Pre-calculate w = exp(-dE/T) for current T
            my_averages[m] = 0.0;               // Initialize local averages to zero.
            averages[m] = 0;                    // Initialize all total averages to zero
        }

        //initialize_spin_config_ordered(spins, num_spins);   // Set initial spin config for current T
        //initialize_spin_config_rng(spins, num_spins);       // Set initial spin config for current T
        compute_energy_and_moment(spins, num_spins, E, M);  // Compute initial energy and moment
        //print_spin_array(spins, num_spins);

        //metropolis(spins, num_spins, num_mc_cycles, delta_E, exp_dE, E, M, expectation_values);
        /* ================================== MC cycles ================================== */
        for (int n = 1; n <= my_num_mc_cycles; n++){    // Loop over number of local MC cycles
            for (int i = 0; i < total_num_spins; i++){        // Loop over rows in spins array
                int k_rng = (int) (uniform(rng)*(double)num_spins);  // Random spin row index
                int l_rng = (int) (uniform(rng)*(double)num_spins);  // Random spin column index
                int dE = 2*spins[k_rng][l_rng]*
                        (spins[k_rng][get_periodic_index(l_rng+1, num_spins)] +
                         spins[get_periodic_index(k_rng+1, num_spins)][l_rng] +
                         spins[k_rng][get_periodic_index(l_rng-1, num_spins)] +
                         spins[get_periodic_index(k_rng-1, num_spins)][l_rng]);
                //std::cout << dE << std::endl;

                /* Metropolis algorithm for accepting or discarding new configuration.
                 * =================================================================== */
                if (dE <= 0){    // If dE <= 0
                    //std::cout << "Flipped a spin" << std::endl;
                    spins[k_rng][l_rng] *= -1;
                    M += (double) (2*spins[k_rng][l_rng]);
                    E += (double) dE;
                    my_accepted_states += 1;
                }

                else {
                    double current_exp_dE = 0.0;    // To hold current w = exp(-dE/T)

                    for (int m = 0; m < 5; m++){    // Loop over delta_E array
                        //std::cout << dE << ", " << delta_E[m] << std::endl;
                        if (dE == delta_E[m]){
                            //std::cout << "Found pre-calculated exp_dE" << std::endl;
                            current_exp_dE = exp_dE[m];     // Found pre-calculated w = exp(-dE/T)
                            break;                          // No need to search anymore if found
                        }
                    }

                    if (uniform(rng) <= current_exp_dE){    // If r <= exp(-dE/T)
                        //std::cout << "Flipped a spin" << std::endl;
                        spins[k_rng][l_rng] *= -1;
                        M += (double) (2*spins[k_rng][l_rng]);
                        E += (double) dE;
                        my_accepted_states += 1;
                    }
                }
                /* =================================================================== */
            }

            /* Update local averages every MC cycle. */
            my_averages[0] += E; my_averages[1] += E*E;
            my_averages[2] += M; my_averages[3] += M*M;
            my_averages[4] += fabs(M);

            //int mc_write = n*num_procs;

//            /* If statement uncommented when counting energy for probability distribution. */
//            if (n >= 6000000 && my_rank == 0){
//                outfile_02 << std::setw(15) << std::setprecision(8) << E/((double)total_num_spins) << std::endl;
//            }

//            /* IF statement uncommented for writing averages as functions of MC cycles. */
//            if (mc_write % 10000 == 0 || n == 1){
//                MPI_Reduce(my_averages, averages, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//                MPI_Reduce(&my_accepted_states, &accepted_states, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

//                if (my_rank == 0){
//                    write_to_file(outfile, T, averages, accepted_states, num_spins, mc_write);  // Write averages to file.
//                }
//            }
        }

        /* =========================== Post MC-cycles logistics ===========================*/
        /* Sum expectation values and num. accepted states from all processes and collect at process 0.*/
        MPI_Reduce(my_averages, averages, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&my_accepted_states, &accepted_states, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (my_rank == 0){
            write_to_file(outfile, T, averages, accepted_states, num_spins, num_mc_cycles);    // Write averages to file.
        }

        T += dT;    // Increase temperature
    }

    /* ============================ Post temperature loop logistics ============================= */
    MPI_Barrier(MPI_COMM_WORLD);    // Sync all processes to get time (below) used by "slowest" process

    if (my_rank == 0){
        time_end = MPI_Wtime();
        time_used = time_end - time_start;
        std::cout << "\nTime used by process " << my_rank << ": " << time_used << " seconds.\n" << std::endl;
        outfile.close();
        outfile_02.close();
    }

    free_array_2D(spins, num_spins);
}


/* Function that initializes 2D spins array with the current spin config except if the
 * temperature is less than 1.5, then it initializes all spins up (+1). */
void initialize_spin_config_ordered(int**& spins, int num_spins){
    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            spins[i][j] = 1;
        }
    }
}


/* Function that initilizes the 2D spin array with random spin values (+1 or -1). */
void initialize_spin_config_rng(int**& spins, int num_spins){
    std::random_device rd;      // Instantiate random device
    std::mt19937_64 rng(rd());  // Instantiate random number generator
    std::uniform_real_distribution<double> uniform(0.0, 1.0);  // Use a uniform dist for the generator

    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            if (uniform(rng) < 0.5){
                spins[i][j] = -1;
            }

            else {
                spins[i][j] = 1;
            }
        }
    }
}

/* Function that computes the energy and magnetic moment of a spin confiuration. */
void compute_energy_and_moment(int** spins, int num_spins, double& E, double& M){
    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            M += (double) spins[i][j];
            E -= (double) (spins[i][j]*(spins[i][get_periodic_index(j+1, num_spins)] +
                    spins[get_periodic_index(i+1, num_spins)][j]));  // Nearest neighbour
        }
    }
}


/* Function that implements periodic boundary conditions, i.e. checks if the proposed_index
 * is out of bounds and, if true, returns the opposite wall index. Otherwise returns input. */
int get_periodic_index(int proposed_index, int array_lenght){
    if (proposed_index < 0){
        return array_lenght - 1;
    }

    else if (proposed_index >= array_lenght){
        return 0;
    }

    else {
        return proposed_index;
    }
}


/* Function that writes num_spins, T, num_mc_cycles, <E>, <|M|>, <M>, <E^2> <M^2>, C_v, X,
 * accepted_states to file (normilized for number of spins). */
void write_to_file(std::ofstream& outfile, double T, double* expectation_values,
                   int accepted_states, int num_spins, int num_mc_cycles){
    double one_over_mcc = 1.0/((double)num_mc_cycles);
    double one_over_total_spins = 1.0/((double)(num_spins*num_spins));

    double E_avg = expectation_values[0]*one_over_mcc;
    double E2_avg = expectation_values[1]*one_over_mcc;
    double M_avg = expectation_values[2]*one_over_mcc;
    double M2_avg = expectation_values[3]*one_over_mcc;
    double Mabs_avg = expectation_values[4]*one_over_mcc;

    double E_var = E2_avg - E_avg*E_avg;
    double M_abs_var = M2_avg - Mabs_avg*Mabs_avg;

    double C_v = E_var*one_over_total_spins/(T*T);
    double X = M_abs_var*one_over_total_spins/T;

    outfile << std::setw(10) << std::setprecision(8) << num_spins;
    outfile << std::setw(18) << std::setprecision(8) << T;
    outfile << std::setw(18) << std::setprecision(8) << num_mc_cycles;
    outfile << std::setw(18) << std::setprecision(8) << E_avg*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << Mabs_avg*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << M_avg*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << E_var*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << M_abs_var*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << C_v;
    outfile << std::setw(18) << std::setprecision(8) << X;
    outfile << std::setw(18) << std::setprecision(8) << accepted_states << std::endl;
}

/* Function that prints out analytical values for the 2x2 lattice case. */
void print_analytical_values(double T){
    double one_over_total_spins = 1.0/(2.0*2.0);
    double k_B = 1.0;
    double beta = 1.0/(k_B*T);

    double Z = 4.0*cosh(8*beta) + 12.0;
    double E = -32.0*sinh(8*beta)/Z;
    double E2 = 256.0*sinh(8*beta)/Z;
    double M = (8.0*exp(8.0*beta) + 16.0)/Z;
    double M2 = 32.0*(exp(8.0*beta) + 1)/Z;
    double C_v = (E2 - E*E)/(k_B*T*T);
    double X = (M2 - M*M)/(k_B*T);

    std::cout << "<E>   = " << E*one_over_total_spins << std::endl;
    std::cout << "<E^2> = " << E2*one_over_total_spins << std::endl;
    std::cout << "<M>   = " << M*one_over_total_spins << std::endl;
    std::cout << "<M^2> = " << M2*one_over_total_spins << std::endl;
    std::cout << "C_v   = " << C_v*one_over_total_spins << std::endl;
    std::cout << "X     = " << X*one_over_total_spins << std::endl << std::endl;
}


/* Function for displaying the 2D spin array to standard out (for debugging). */
void print_spin_array(int** spins, int num_spins){
    std::cout << std::endl;

    for(int i = 0; i < num_spins; i++){
        for(int j = 0; j < num_spins; j++){
            std::cout << std::setw(2) << spins[i][j] << " ";
        }

        std::cout  << std::endl;
    }
}
