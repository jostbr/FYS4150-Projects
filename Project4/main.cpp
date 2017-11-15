
# include <iostream>
# include "ising.h"

int main(){
    double T_0 = 1.0;
    double T_N = 3.0;
    double dT = 0.05;
    int num_spins = 2;
    int num_mc_cycles = 1000000;
    ising(T_0, T_N, dT, num_spins, num_mc_cycles);
}
