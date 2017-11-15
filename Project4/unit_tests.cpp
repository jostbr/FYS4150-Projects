
# include <unit_tests.h>

/* Unit test implementation of periodic boundary condition function. Reports
 * if the lower or upper boundary condition fails, or if both tests are passed. */
void TEST_get_index(){
    bool test_01 = get_index(-1, 10) == 9;
    bool test_02 = get_index(10, 10) == 0;

    if (!test_01 || !test_02){
        if (!test_01){
            std::cout << "TEST FAILED: Lower periodic boundary conditions" << std::endl;
        }

        if (!test_02){
            std::cout << "TEST FAILED: Upper periodic boundary conditions" << std::endl;
        }
    }

    else {
        std::cout << "TEST PASSED: Periodic boundary conditions" << std::endl;
    }
}
