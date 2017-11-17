#ifndef ISING_H
#define ISING_H

# include <iostream>
# include <iomanip>
# include <fstream>
# include <random>
# include <cmath>
# include <ctime>

void ising(double T_min, double T_max, double dT, int num_spins, int num_mc_cycles);
void ising_mpi(double T_0, double T_N, double dT, int num_spins, int num_mc_cycle, int my_rank, int num_procs);
void metropolis(int**& spins, int num_spins, int num_mc_cycles, double* delta_E,
                double* exp_dE, double& E, double& M, double* expectation_values);
void initialize_spin_config(int**& spins, int num_spins, double& E, double& M, double T);
void initialize_spin_config_prev(int**& spins, int num_spins, double T);
void initialize_spin_config_rng(int**& spins, int num_spins);
void compute_energy_and_moment(int** spins, int num_spins, double& E, double& M);
int get_periodic_index(int proposed_index, int array_lenght);
void write_to_file(std::ofstream& outfile, double T, double* expectation_values,
                   int num_spins, int num_mc_cycles);
void print_analytical_values(double T);
void print_spin_array(int** spins, int num_spins);

#endif // ISING_H
