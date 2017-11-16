
# include "ising.h"

/* Function performing the Monte Carlo simulation in the Ising model. Takes input initial and final
 * temperature of the system, a temperature step, number of spins in one direction (of the 2D square
 * lattice) and the number of Monte Carlo cycle one wishes to perform (the more the better). */
void ising(double T_min, double T_max, double dT, int num_spins, int num_mc_cycles){
    double E, M;                    // To hold energy and magentic moment each temperature iteration
    double T = T_min;               // Temperature starts at minimum
    double delta_E[5], exp_dE[5];   // Arrays to hold pre-caclualted energy differences
    double expectation_values[5];   // Array to hold expectation values of E, E^2, M, M^2, |M|

    std::random_device rd;      // Instantiate random device
    std::mt19937_64 rng(rd());  // Instantiate random number generator
    std::uniform_real_distribution<double> uniform(0.0, 1.0);  // Use a uniform dist for the generator

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for(int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    delta_E[0] = -8.0; delta_E[1] = -4.0; delta_E[2] = -0.0; delta_E[3] = 4.0; delta_E[4] = 8.0;  // Possible dE's

    /* Instantiate output file object and write header info to file. */
    std::ofstream outfile;
    outfile.open("results.txt");
    outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    outfile << std::setw(18) << "Num spins" << std::setw(18) << "MC cycles" << std::setw(18)
            << "Temperature" << std::setw(18) << "Mean energy" << std::setw(18) << "Mean abs mag"
            << std::setw(18) << "Heat cap" << std::setw(18) << "Mag susc" << std::endl;

    clock_t t_0 = clock();

    /* Main loop over temperature doing Monte Carlo simulations each temperature step. */
    while (T <= T_max){
        E = 0;      // Start energy at zero for every temperature
        M = 0;      // Start magnetic moment at zero for every temperature

        for (int m = 0; m < 5; m++){
            exp_dE[m] = exp(-delta_E[m]/T);     // Initialize w = exp(-dE/T) for current T
            expectation_values[m] = 0.0;        // Initialize all averages to zero
        }

        //initialize_spin_config_prev(spins, num_spins, T);  // Set initial spin config for current T
        initialize_spin_config_rng(spins, num_spins);  // Set initial spin config for current T
        compute_energy_and_moment(spins, num_spins, E, M);
        print_spin_array(spins, num_spins);

        //metropolis(spins, num_spins, num_mc_cycles, delta_E, exp_dE, E, M, expectation_values);
        for (int n = 0; n < num_mc_cycles; n++){    // Loop over number of MC cycles
            for (int i = 0; i < num_spins; i++){        // Loop over rows in spins array
                for (int j = 0; j < num_spins; j++){        // Loop over columns in spins array
                    int k_rng = (int) uniform(rng)*(double)num_spins;  // Random spin row index
                    int l_rng = (int) uniform(rng)*(double)num_spins;  // Random spin column index
                    int dE = (double) 2*spins[k_rng][l_rng]*
                            (spins[k_rng][get_periodic_index(l_rng+1, num_spins)] +
                             spins[get_periodic_index(k_rng+1, num_spins)][l_rng] +
                             spins[k_rng][get_periodic_index(l_rng-1, num_spins)] +
                             spins[get_periodic_index(k_rng-1, num_spins)][l_rng]);



                    if (dE <= 0){    // If dE <= 0 or if random(0,1) < exp(-dE/T)
                        //std::cout << "Flipped a spin" << std::endl;
                        spins[k_rng][l_rng] *= -1;
                        M += (double) 2*spins[k_rng][l_rng];
                        E += (double) dE;
                    }

                    else {
                        double current_exp_dE = 0.0;    // To hold current w = exp(-dE/T)

                        for (int m = 0; m < 5; m++){
                            if ((double) dE == delta_E[m]){
                                //std::cout << "Found pre-calculated exp_dE" << std::endl;
                                current_exp_dE = exp_dE[m];     // Found pre-calculated w = exp(-dE/T)
                            }
                        }

                        if (uniform(rng) <= current_exp_dE){
                            spins[k_rng][l_rng] *= -1;
                            M += (double) 2*spins[k_rng][l_rng];
                            E += (double) dE;
                        }
                    }
                }
            }

            expectation_values[0] += E; expectation_values[1] += E*E;
            expectation_values[2] += M; expectation_values[3] += M*M;
            expectation_values[4] += fabs(M);
        }

        write_to_file(outfile, T, expectation_values, num_spins, num_mc_cycles);    // Write averages to file.

        T += dT;
    }

    clock_t t_1 = clock();
    double time_used = (double)(t_1 - t_0)/CLOCKS_PER_SEC;
    std::cout << "\nTime used: " << std::setprecision(8) << time_used << " seconds.\n" << std::endl;

    outfile.close();

    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++){
        delete [] spins[i];
    }

    delete [] spins;
}

//void metropolis(int**& spins, int num_spins, int num_mc_cycles, double* delta_E, double* exp_dE, double& E, double& M, double* expectation_values){
//    std::random_device rd;      // Instantiate random device
//    std::mt19937_64 rng(rd());  // Instantiate random number generator
//    std::uniform_real_distribution<double> uniform(0.0, 1.0);  // Use a uniform dist for the generator

//    for (int n = 0; n < num_mc_cycles; n++){    // Loop over number of MC cycles
//        for (int i = 0; i < num_spins; i++){        // Loop over rows in spins array
//            for (int j = 0; j < num_spins; j++){        // Loop over columns in spins array
//                int k_rng = (int) uniform(rng)*(double)num_spins;  // Random spin row index
//                int l_rng = (int) uniform(rng)*(double)num_spins;  // Random spin column index
//                int dE = (double) 2*spins[k_rng][l_rng]*
//                        (spins[k_rng][get_periodic_index(l_rng+1, num_spins)] +
//                         spins[get_periodic_index(k_rng+1, num_spins)][l_rng] +
//                         spins[k_rng][get_periodic_index(l_rng-1, num_spins)] +
//                         spins[get_periodic_index(k_rng-1, num_spins)][l_rng]);

//                double current_exp_dE = 0.0;    // To hold current w = exp(-dE/T)

//                for (int m = 0; m < 5; m++){
//                    if ((double) dE == delta_E[m]){
//                        //std::cout << "Found pre-calculated exp_dE" << std::endl;
//                        current_exp_dE = exp_dE[m];     // Found pre-calculated w = exp(-dE/T)
//                    }
//                }

//                if (dE <= 0 || uniform(rng) <= current_exp_dE){    // If dE <= 0 or if random(0,1) < exp(-dE/T)
//                    //std::cout << "Flipped a spin" << std::endl;
//                    spins[k_rng][l_rng] *= -1;
//                    M += (double) 2*spins[k_rng][l_rng];
//                    E += (double) dE;
//                }
//            }
//        }

//        expectation_values[0] += E; expectation_values[1] += E*E;
//        expectation_values[2] += M; expectation_values[3] += M*M;
//        expectation_values[4] += fabs(M);
//    }
//}


/* Function that initializes 2D spins array with the current spin config except if the
 * temperature is less than 1.5, then it initializes all spins up (+1). */
void initialize_spin_config_prev(int**& spins, int num_spins, double T){
    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            if (T < 1.5){
                spins[i][j] = 1;
            }
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
            M += (double)spins[i][j];
            E -= (double)spins[i][j]*(spins[i][get_periodic_index(j+1, num_spins)] +
                    spins[get_periodic_index(i+1, num_spins)][j]);  // Nearest neighbour
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


/* Function that writes T, E, E^2, M M^2 and |M| to file (normilized for number of spins). */
void write_to_file(std::ofstream& outfile, double T, double* expectation_values, int num_spins, int num_mc_cycles){
    double one_over_mcs = 1.0/((double)num_mc_cycles);
    double one_over_total_spins = 1.0/((double)(num_spins*num_spins));

    double E_avg = expectation_values[0]*one_over_mcs;
    double E2_avg = expectation_values[1]*one_over_mcs;
    double M_avg = expectation_values[2]*one_over_mcs;
    double M2_avg = expectation_values[3]*one_over_mcs;
    double Mabs_avg = expectation_values[4]*one_over_mcs;

    double C_v = (E2_avg - E_avg*E_avg)*one_over_total_spins/(T*T);
    double X = (M2_avg - Mabs_avg*Mabs_avg)*one_over_total_spins/T;

    outfile << std::setw(18) << std::setprecision(8) << num_spins;
    outfile << std::setw(18) << std::setprecision(8) << num_mc_cycles;
    outfile << std::setw(18) << std::setprecision(8) << T;
    outfile << std::setw(18) << std::setprecision(8) << E_avg*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << Mabs_avg*one_over_total_spins;
    outfile << std::setw(18) << std::setprecision(8) << C_v;
    outfile << std::setw(18) << std::setprecision(8) << X << std::endl;
}

/* Function that prints out analytical values for the 2x2 lattice case. */
void print_analytical_values(int num_spins, double T){
    double one_over_total_spins = 1.0/((double)num_spins*num_spins);
    double k_B = 1.0;
    double beta = 1/(k_B*T);

    double Z = 4*cosh(8*beta) + 12;
    double E = -32*sinh(8*beta)/Z;
    double E2 = -256*sinh(8*beta)/Z;
    double M = (8*exp(8*beta) + 16)/Z;
    double M2 = 32*(32*exp(8*beta) + 32)/Z;
    double C_v = (E2 - E*E)/(k_B*T*T);
    double X = (M2 - M*M)/(k_B*T);

    std::cout << "<E> = " << E*one_over_total_spins << std::endl;
    std::cout << "<E^2> = " << E2*one_over_total_spins << std::endl;
    std::cout << "<M> = " << M*one_over_total_spins << std::endl;
    std::cout << "<M^2> = " << M2*one_over_total_spins << std::endl;
    std::cout << "C_v = " << C_v*one_over_total_spins << std::endl;
    std::cout << "X = " << X*one_over_total_spins << std::endl;
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
