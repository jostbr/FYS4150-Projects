
# include <iostream>
# include "periodic_solver_1d.hpp"
# include "basin_solver_1d.hpp"
# include "array_alloc.hpp"
# include "unit_tests.hpp"

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag_general();
    TEST_set_basin_IC();
    std::cout << "===============================================\n" << std::endl;

    int N = 100;
    double* init_con;
    alloc_array_1D(init_con, N);

    for (int i = 0; i < N; i++){
        init_con[i] = 0.0;
    }

    basin_solver_1d solver(4.0, 1.2, N, 500);
    solver.set_boundary_conditions(0.0, 0.0);
    solver.set_initial_condition(init_con);

    periodic_solver_1d psolver(1.0,1.0,10,100);
    init_con[0] = 700.0;
    psolver.set_initial_condition(init_con);



    return 0;
}
