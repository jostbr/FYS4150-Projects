
# include <iostream>
# include "periodic_solver_1d.hpp"
# include "basin_solver_1d.hpp"
# include "unit_tests.hpp"

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag();
    std::cout << "===============================================" << std::endl;


    basin_solver_1d solver(4.0, 1.2, 100, 500);
    //std::cout << solver.T << std::endl;



    return 0;
}
