
# include <iostream>
# include <armadillo>
# include "unit_tests.h"
# include "jacobi.h"

int main(){
    std::cout << "EXECUTING UNIT TESTS..." << std::endl
              << "============================================" << std::endl;
    test_max_non_diag();
    // Do some more tests
    std::cout << "============================================" << std::endl;

    return 0;
}

