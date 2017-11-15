#ifndef ISING_H
#define ISING_H

# include <iostream>
# include <fstream>
# include <random>
# include <cmath>
# include <ctime>

void ising(double T_0, double T_N, double dT, int num_spins, int num_mc_cycle);
void initialize_spin_config(int**& spins, int num_spins, double& E, double& M, double T);
int get_periodic_index(int proposed_index, int array_lenght);
void print_spin_array(int** spins, int num_spins);

#endif // ISING_H
