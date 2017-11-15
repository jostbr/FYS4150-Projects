
# include "ising.h"

/* Function performing the Monte Carlo simulation in the Ising model. Takes input initial and final
 * temperature of the system, a temperature step, number of spins in one direction (of the 2D square
 * lattice) and the number of Monte Carlo cycle one wishes to perform (the more the better). */
void ising(double T_min, double T_max, double dT, int num_spins, int num_mc_cycle){
    double E, M;                    // To hold energy and magentic moment each temperature iteration
    double T = T_min;               // Temperature starts at minimum
    double delta_E[5];              // Static array to hold pre-caclualted energy differences
    double expectation_values[5];   // Array to hold expectation values of E, E^2, M, M^2, |M|

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for(int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    print_spin_array(spins, num_spins);
    //initialize_config(spins, num_spins, E, M, T);
    print_spin_array(spins, num_spins);

    /* Main loop over temperature doing Monte Carlo simulations each step. */
    while (T <= T_max){
        E = 0;
        M = 0;

        T += dT;
    }

    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++) {
        delete [] spins[i];
    }

    delete [] spins;
}

void initialize_config(int**& spins, int num_spins, double& E, double& M, double T){
    for (int i = 0; i < num_spins; i++){
        for (int j = 0; i < num_spins; j++){
            if (T < 1.5){
                spins[i][j] = 1;
            }

            M += (double)spins[i][j];
        }
    }
}

/* Function that implements periodic boundary conditions, i.e. checks if the proposed_index
 * is out of bounds and, if true, returns the opposite wall index. Otherwise returns input. */
int get_index(int proposed_index, int num_spins){
    if (proposed_index < 0){
        return num_spins - 1;
    }

    else if (proposed_index >= num_spins){
        return 0;
    }

    else {
        return proposed_index;
    }
}

/* Function for displaying the 2D spin array to standard out (for debugging). */
void print_spin_array(int** spins, int num_spins){
    for(int i = 0; i < num_spins; i++){
        for(int j = 0; j < num_spins; j++){
            std::cout << spins[i][j] << " ";
        }

        std::cout  << std::endl;
    }
}
