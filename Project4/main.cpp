
# include <iostream>
# include "ising.h"
# include "unit_tests.h"

int main(){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_get_index();
    std::cout << "===============================================" << std::endl;

    double T_0 = 1.0;
    double T_N = 3.0;
    double dT = 0.05;
    int num_spins = 2;
    int num_mc_cycles = 1000000;
    ising(T_0, T_N, dT, num_spins, num_mc_cycles);
}
