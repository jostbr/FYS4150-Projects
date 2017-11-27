
# include <iostream>
# include <string>
# include <cmath>
# include "periodic_solver_1d.hpp"
# include "basin_solver_1d.hpp"
# include "array_alloc.hpp"
# include "unit_tests.hpp"

void func(double *arr, int n);

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag_general();
    TEST_tridiag_ferrari();
    TEST_set_initial_condition_basin();
    TEST_set_initial_condition_periodic();
    std::cout << "===============================================\n" << std::endl;


    double pi = acos(-1.0);

    double L = 1.0;         // Length of domain (non-dim)
    double dx = 1.0/40.0;   // Spatial step
    int N = L/dx + 1;       // Number of spatial points

    double T = 100.0;        // Upper time limit
    double dt = 0.001;     // Temporal step

    double *init_psi, *init_zeta;   // Initial conditions
    alloc_array_1D(init_psi, N);
    alloc_array_1D(init_zeta, N);

    double x;

    for (int i = 0; i < N; i++){
        x = i*dx;
        //std::cout << x << std::endl;
        init_psi[i] = sin(4.0*pi*x);
        init_zeta[i] = -16.0*pi*pi*sin(4.0*pi*x);
    }

    std::string fileout_01 = "sine_basin_euler.txt";
    basin_solver_1d be_solver(dx, N, dt, T, fileout_01);
    be_solver.set_initial_condition(init_psi, init_zeta);
    be_solver.basin_euler();

    std::string fileout_02 = "sine_basin_leapfrog.txt";
    basin_solver_1d bl_solver(dx, N, dt, T, fileout_02);
    bl_solver.set_initial_condition(init_psi, init_zeta);
    bl_solver.basin_leapfrog();

    std::string fileout_03 = "sine_periodic_euler.txt";
    periodic_solver_1d pe_solver(dx, N, dt, T, fileout_03);
    pe_solver.set_initial_condition(init_psi, init_zeta);
    pe_solver.periodic_euler();

    std::string fileout_04 = "sine_periodic_leapfrog.txt";
    periodic_solver_1d pb_solver(dx, N, dt, T, fileout_04);
    pb_solver.set_initial_condition(init_psi, init_zeta);
    pb_solver.periodic_leapfrog();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);

    return 0;
}
