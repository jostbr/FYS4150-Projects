
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
    //TEST_set_basin_IC();
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

    std::string fileout = "results.txt";
    basin_solver_1d solver(dx, dt, N, T, fileout);
    solver.set_boundary_conditions(0.0, 0.0);
    solver.set_initial_condition(init_psi, init_zeta);
    solver.basin_euler();

    free_array_1D(init_psi);
    free_array_1D(init_zeta);

    double* arr;
    alloc_array_1D(arr, 6);

    for (int i = 0; i < 6; i++){
        arr[i] = i;
    }

    func(arr, 6);

    for (int i = 0; i < 6; i++){
        std::cout << "\nGadhjgsfv" << std::endl;
        std::cout << arr[i] << std::endl;
    }

    free_array_1D(arr);

    return 0;
}

void func(double* arr, int n){
    for (int i = 0; i < n; i++){
        std::cout << "Gadhjgsfv" << std::endl;
        std::cout << arr[i] << std::endl;
    }

    arr[0] = 276.0;
}
