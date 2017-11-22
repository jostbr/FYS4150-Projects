
# include <iostream>
# include "unit_tests.hpp"

int main(int argc, char* argv[]){
    std::cout << "\nRUNNING TESTS..." << std::endl;
    std::cout << "===============================================" << std::endl;
    TEST_tridiag();
    std::cout << "===============================================" << std::endl;

    return 0;
}
