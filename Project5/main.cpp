
# include <iostream>
# include <string>
# include <cmath>
# include "basin_solver_1d.hpp"
# include "periodic_solver_1d.hpp"
# include "basin_solver_2d.hpp"
# include "periodic_solver_2d.hpp"
# include "array_alloc.hpp"
# include "unit_tests.hpp"

void run_1d_simulations();
void run_2d_simulations();
void func(double *arr, int n);

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag_general();
    TEST_tridiag_ferrari();
    TEST_set_initial_condition_basin();
    TEST_set_initial_condition_periodic();
    std::cout << "===============================================\n" << std::endl;

    //run_1d_simulations();
    run_2d_simulations();

    std::cout << "Sup?" << std::endl;

    return 0;
}

void run_1d_simulations(){
    /* Define some simulation parameters. */
    double pi = acos(-1.0);

    double L = 1.0;         // Length of domain (non-dim)
    double dx = 1.0/40.0;   // Spatial step in x-direction
    int N_x = L/dx + 1;       // Number of spatial points in x-direction

    double T = 50.0;        // Upper time limit
    double dt = 0.001;      // Time step

    /* Initial conditions for 1D case. */
    double *init_psi, *init_zeta;   // Initial conditions
    alloc_array_1D(init_psi, N_x);
    alloc_array_1D(init_zeta, N_x);

    double x;

    for (int i = 0; i < N_x; i++){
        x = i*dx;
        //std::cout << x << std::endl;
        init_psi[i] = sin(4.0*pi*x);
        init_zeta[i] = -16.0*pi*pi*sin(4.0*pi*x);
    }

//    std::string fileout_01 = "sine_basin_euler.txt";
//    basin_solver_1d be_solver(dx, N_x, dt, T, fileout_01);
//    be_solver.set_initial_condition(init_psi, init_zeta);
//    be_solver.basin_euler();

//    std::string fileout_02 = "sine_basin_leapfrog.txt";
//    basin_solver_1d bl_solver(dx, N_x, dt, T, fileout_02);
//    bl_solver.set_initial_condition(init_psi, init_zeta);
//    bl_solver.basin_leapfrog();

    std::string fileout_03 = "sine_periodic_euler.txt";
    periodic_solver_1d pe_solver(dx, N_x, dt, T, fileout_03);
    pe_solver.set_initial_condition(init_psi, init_zeta);
    pe_solver.periodic_euler();

//    std::string fileout_04 = "sine_periodic_leapfrog.txt";
//    periodic_solver_1d pb_solver(dx, N_x, dt, T, fileout_04);
//    pb_solver.set_initial_condition(init_psi, init_zeta);
//    pb_solver.periodic_leapfrog();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);
}

void run_2d_simulations(){
    /* Define some simulation parameters. */
    double pi = acos(-1.0);

    double L = 1.0;         // Length of domain (non-dim)
    double dx = 1.0/40.0;   // Spatial step in x-direction
    double dy = 1.0/40.0;   // Spatial step in x-direction
    int N_x = L/dx + 1;       // Number of spatial points in x-direction
    int N_y = L/dy + 1;       // Number of spatial points in x-direction

    double T = 100.0;        // Upper time limit
    double dt = 0.001;      // Time step

    /* Initial conditions. */
    double *init_psi, *init_zeta;   // Initial conditions
    alloc_array_1D(init_psi, N_x*N_y);
    alloc_array_1D(init_zeta, N_x*N_y);

    double x, y;

    for (int i = 0; i < N_x; i++){
        for (int j = 0; j < N_x; j++){
            x = i*dx;
            y = j*dy;
            init_psi[i*N_y + j] = sin(pi*y)*sin(4.0*pi*x);
            init_zeta[i*N_y + j] = -17.0*pi*pi*sin(pi*y)*sin(4.0*pi*x);
        }
    }

//    std::string fileout_05 = "results_2d.txt";
//    basin_solver_2d b2d_solver(dx, dy, N_x, N_y, dt, T, fileout_05);
//    b2d_solver.set_initial_condition(init_psi, init_zeta);
//    b2d_solver.basin_leapfrog();

    std::string fileout_06 = "results_2d.txt";
    periodic_solver_2d p2d_solver(dx, dy, N_x, N_y, dt, T, fileout_06);
    p2d_solver.set_initial_condition(init_psi, init_zeta);
    p2d_solver.periodic_leapfrog();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);
}

