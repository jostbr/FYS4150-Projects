
# include <iostream>
# include <mpi.h>
# include "ising.h"
# include "ising_mpi.h"
# include "unit_tests.h"

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_get_index();
    TEST_compute_energy_and_moment();
    TEST_initialize_spin_config_prev();
    TEST_initialize_spin_config_rng();
    std::cout << "===============================================" << std::endl;

    std::cout << "\nAnalytical values for 2x2 lattice (per spin) with T = 1.0" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    print_analytical_values(1.0);

//    double T_0 = 1.0;
//    double T_N = 1.0;
//    double dT = 0.05;
//    int num_spins = 2;
//    int num_mc_cycles = 1000000;
//    ising(T_0, T_N, dT, num_spins, num_mc_cycles);
    ising_mpi();

    int my_rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    std::cout << "Hello world from process " << my_rank << std::endl;
    MPI_Finalize();
}
