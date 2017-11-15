
# include "ising.h"

/* Function performing the Monte Carlo simulation in the Ising model. Takes input initial and final
 * temperature of the system, a temperature step, number of spins in one direction (of the 2D square
 * lattice) and the number of Monte Carlo cycle one wishes to perform (the more the better). */
void ising(double T_min, double T_max, double dT, int num_spins, int num_mc_cycles){
    double E, M;                    // To hold energy and magentic moment each temperature iteration
    double T = T_min;               // Temperature starts at minimum
    double delta_E[5], exp_dE[5];   // Arrays to hold pre-caclualted energy differences
    double expectation_values[5];   // Array to hold expectation values of E, E^2, M, M^2, |M|

    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for(int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    delta_E[0] = -8.0; delta_E[1] = -4.0; delta_E[2] = -0.0; delta_E[3] = 4.0; delta_E[4] = 8.0;



    /* Main loop over temperature doing Monte Carlo simulations each step. */
    while (T <= T_max){
        E = 0;      // Start energy at zero for every
        M = 0;

        for (int m = 0; m < 5; m++){
            exp_dE[m] = exp(delta_E[m]/T);  // Initialize w = exp(dE/T) for current T
            expectation_values[m] = 0.0;    // Initialize all averages to zero
        }

        initialize_spin_config(spins, num_spins, E, M, T);  // Set initial spin config for current T
        //print_spin_array(spins, num_spins);

        for (int n = 0; n < num_mc_cycles; n++){
            for (int i = 0; i < num_spins; i++){
                for (int j = 0; j < num_spins; j++){
                    int k_rng = (int) dist(rng)*(double)num_spins;
                    int l_rng = (int) dist(rng)*(double)num_spins;
                    int dE = (double) spins[k_rng][l_rng]*
                            (spins[k_rng][get_periodic_index(l_rng+1, num_spins)] +
                             spins[get_periodic_index(k_rng+1, num_spins)][l_rng] +
                             spins[k_rng][get_periodic_index(l_rng-1, num_spins)] +
                             spins[get_periodic_index(k_rng-1, num_spins)][l_rng]); // Need to explain this formula

                    double current_exp_dE = 0.0;

                    for (int m = 0; m < 5; m++){
                        if (dE == exp_dE[m]){
                            std::cout << "Found pre-calculated exp_dE" << std::endl;
                            current_exp_dE = dE;
                        }
                    }

                    if (delta_E <= 0 || dist(rng) <= current_exp_dE){
                        std::cout << "Found pre-calculated exp_dE" << std::endl;
                        spins[k_rng][l_rng] *= -1;
                        M += (double) 2*spins[k_rng][l_rng];
                        E += (double) dE;
                    }
                }
            }

            expectation_values[0] += E; expectation_values[1] += E*E;
            expectation_values[2] += M; expectation_values[3] += M*M;
            expectation_values[4] += fabs(M);
        }

        T += dT;
    }

    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++){
        delete [] spins[i];
    }

    delete [] spins;
}

/* Function that initializes 2D spins array with a spin configuration and computes the
 * corresponding energy and magnetic moment of the configuration. */
void initialize_spin_config(int**& spins, int num_spins, double& E, double& M, double T){
    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            if (T < 1.5){
                spins[i][j] = 1;
            }
        }
    }

    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            M += (double)spins[i][j];
            E -= (double)spins[i][j]*(spins[i][get_periodic_index(j+1, num_spins)] +
                    spins[get_periodic_index(i+1, num_spins)][j]);
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

/* Function for displaying the 2D spin array to standard out (for debugging). */
void print_spin_array(int** spins, int num_spins){
    std::cout << std::endl;

    for(int i = 0; i < num_spins; i++){
        for(int j = 0; j < num_spins; j++){
            std::cout << spins[i][j] << " ";
        }

        std::cout  << std::endl;
    }
}
