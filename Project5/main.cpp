
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

int main(){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag_general();
    TEST_tridiag_ferrari();
    TEST_set_initial_condition_basin();
    TEST_set_initial_condition_periodic();
    std::cout << "===============================================\n" << std::endl;

    run_1d_simulations();
    //run_2d_simulations();

    std::cout << "Sup?" << std::endl;

    return 0;
}

void run_1d_simulations(){
    /* Define some simulation parameters. */
    double pi = acos(-1.0);

    double L = 1.0;         // Length of domain (non-dim)
    double dx = 1.0/40.0;   // Spatial step in x-direction
    int N_x = L/dx + 1;       // Number of spatial points in x-direction

    double T = 200.0;        // Upper time limit
    double dt = 0.025;      // Time step

    double sigma = 0.1;     // For Gaussian initial condition

    /* Initial conditions for 1D case. */
    double *init_psi, *init_zeta, *init_psi_gaussian, *init_zeta_gaussian;   // Initial conditions
    alloc_array_1D(init_psi, N_x);
    alloc_array_1D(init_zeta, N_x);
    alloc_array_1D(init_psi_gaussian, N_x);
    alloc_array_1D(init_zeta_gaussian, N_x);

    double x;

    for (int i = 0; i < N_x; i++){
        x = i*dx;
        //std::cout << x << std::endl;
        init_psi[i] = sin(4.0*pi*x);
        init_zeta[i] = -16.0*pi*pi*sin(4.0*pi*x);

        init_psi_gaussian[i] = exp(-((x-0.5)/sigma)*((x-0.5)/sigma));
        init_zeta_gaussian[i] = (4.0*((x - 0.5)/(sigma*sigma))*((x - 0.5)/(sigma*sigma)) - 2.0/(sigma*sigma))
                *exp(-((x-0.5)/sigma)*((x-0.5)/sigma));
    }

    /* Comment or decomment blocks below to choose what case to run. */
    /* ============================================================= */

//    std::string fileout_01 = "sine_basin_euler.txt";
//    basin_solver_1d be_solver(dx, N_x, dt, T, fileout_01);
//    be_solver.set_initial_condition(init_psi, init_zeta);
//    be_solver.basin_euler();

    std::string fileout_02 = "sine_basin_leapfrog.txt";
    basin_solver_1d bl_solver(dx, N_x, dt, T, fileout_02);
    bl_solver.set_initial_condition(init_psi, init_zeta);
    bl_solver.basin_leapfrog();

//    std::string fileout_03 = "sine_periodic_euler.txt";
//    periodic_solver_1d pe_solver(dx, N_x, dt, T, fileout_03);
//    pe_solver.set_initial_condition(init_psi, init_zeta);
//    pe_solver.periodic_euler();

//    std::string fileout_04 = "sine_periodic_leapfrog.txt";
//    periodic_solver_1d pl_solver(dx, N_x, dt, T, fileout_04);
//    pl_solver.set_initial_condition(init_psi, init_zeta);
//    pl_solver.periodic_leapfrog();

//    std::string fileout_06 = "gauss_periodic_leapfrog.txt";
//    periodic_solver_1d pl_solver_2(dx, N_x, dt, T, fileout_06);
//    pl_solver_2.set_initial_condition(init_psi_gaussian, init_zeta_gaussian);
//    pl_solver_2.periodic_leapfrog();

//    init_psi_gaussian[0] = 0.0;
//    init_psi_gaussian[N_x-1] = 0.0;
//    std::string fileout_05 = "gauss_basin_leapfrog.txt";
//    basin_solver_1d bl_solver_2(dx, N_x, dt, T, fileout_05);
//    bl_solver_2.set_initial_condition(init_psi_gaussian, init_zeta_gaussian);
//    bl_solver_2.basin_leapfrog();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);
    free_array_1D(init_psi_gaussian);
    free_array_1D(init_zeta_gaussian);
}

void run_2d_simulations(){
    /* Define some simulation parameters. */
    double pi = acos(-1.0);

    double L = 1.0;         // Length of domain (non-dim)
    double dx = 1.0/40.0;   // Spatial step in x-direction
    double dy = 1.0/40.0;   // Spatial step in x-direction
    int N_x = L/dx + 1;       // Number of spatial points in x-direction
    int N_y = L/dy + 1;       // Number of spatial points in x-direction

    double T = 200.0;        // Upper time limit
    double dt = 0.01;      // Time step

    /* Initial conditions. */
    double *init_psi, *init_zeta;   // Initial conditions
    alloc_array_1D(init_psi, N_x*N_y);
    alloc_array_1D(init_zeta, N_x*N_y);

    double x, y;

    for (int i = 0; i < N_x; i++){
        for (int j = 0; j < N_x; j++){
            x = i*dx;
            y = j*dy;
            init_psi[i*N_y + j] = sin(4.0*pi*x);
            init_zeta[i*N_y + j] = -16.0*pi*pi*sin(4.0*pi*x);
        }
    }

    /* Comment or decomment blocks below to choose what case to run. */
    /* ============================================================= */

//    std::string fileout_05 = "results_2d.txt";
//    basin_solver_2d b2d_solver(dx, dy, N_x, N_y, dt, T, fileout_05);
//    b2d_solver.set_initial_condition(init_psi, init_zeta);
//    b2d_solver.basin_leapfrog();

    std::string fileout_06 = "periodic_sine_2d.txt";
    periodic_solver_2d p2d_solver(dx, dy, N_x, N_y, dt, T, fileout_06);
    p2d_solver.set_initial_condition(init_psi, init_zeta);
    p2d_solver.periodic_leapfrog();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);
}

