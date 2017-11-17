
# include <unit_tests.h>

/* Unit test for implementation of periodic boundary condition function. Reports
 * if the lower or upper boundary condition fails, or if both tests are passed. */
void TEST_get_index(){
    bool test_01 = get_periodic_index(-1, 10) == 9;
    bool test_02 = get_periodic_index(10, 10) == 0;

    if (!test_01 || !test_02){
        if (!test_01) std::cout << "TEST FAILED: Lower periodic boundary conditions" << std::endl;
        if (!test_02) std::cout << "TEST FAILED: Upper periodic boundary conditions" << std::endl;
    }

    else std::cout << "TEST PASSED: Periodic boundary conditions" << std::endl;
}

/* Unit test for spin configuration initialization function. Uses a known 2x2 lattice (all spin +1)and
 * checks if the function initializes the spins correctly and produces the expected energy and mag. moments. */
void TEST_compute_energy_and_moment(){
    double E = 0.0, M = 0.0;    // t hold the computed energy and magnetic moments of the config
    double T = 1.0;             // Sample temperature (below 1.5 to make sure all spins are +1)
    int num_spins = 2;          // Size of lattice, i.e. 2x2

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for (int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    initialize_spin_config_prev(spins, num_spins, T);   // Call spin config initializer (all spins +1)
    compute_energy_and_moment(spins, num_spins, E, M);  // Compute energy and moment for the config

    /* Compare with analytical values for E and M for the 2x2 lattice. */
    if (E != -8.0 || M != 4.0){
        if (E != -8.0) std::cout << "TEST FAILED: Initialize energy" << std::endl;
        if (M != 4.0) std::cout << "TEST FAILED: Initialize magnetic moment" << std::endl;
    }

    else std::cout << "TEST PASSED: Energy and magnetic moment" << std::endl;


    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++){
        delete [] spins[i];
    }

    delete [] spins;
}

/* Unit test that checks if the random spins config initializer gives spins +1 or -1 for all spins. */
void TEST_initialize_spin_config_prev(){
    double T = 1.0;
    int num_spins = 2;          // Size of lattice, i.e. 2x2

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for (int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    initialize_spin_config_prev(spins, num_spins, T);   // Initialize all +1 spin config

    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            if (spins[i][j] != 1){
                std::cout << "TEST FAILED: All +1 spin configuration" << std::endl;
                return;
            }
        }
    }

    std::cout << "TEST PASSED: All +1 spin configuration" << std::endl;

    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++){
        delete [] spins[i];
    }

    delete [] spins;
}

/* Unit test that checks if the random spins config initializer gives spins +1 or -1 for all spins. */
void TEST_initialize_spin_config_rng(){
    int num_spins = 2;          // Size of lattice, i.e. 2x2

    /* Allocate memory for 2D spin array. */
    int** spins = new int*[num_spins];
    for (int i = 0; i < num_spins; i++){
        spins[i] = new int[num_spins];
    }

    initialize_spin_config_rng(spins, num_spins);   // Initialize random spin config

    for (int i = 0; i < num_spins; i++){
        for (int j = 0; j < num_spins; j++){
            if (spins[i][j] != -1 && spins[i][j] != 1){
                std::cout << "TEST FAILED: Random spin configuration" << std::endl;
                return;
            }
        }
    }

    std::cout << "TEST PASSED: Random spin configuration" << std::endl;

    /* Deallocate memory for spin array. */
    for(int i = 0; i < num_spins; i++){
        delete [] spins[i];
    }

    delete [] spins;
}

