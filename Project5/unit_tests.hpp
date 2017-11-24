
#ifndef UNIT_TESTS_HPP
#define UNIT_TESTS_HPP

# include <iostream>
# include <string>
# include "poisson.hpp"
# include "basin_solver_1d.hpp"
# include "periodic_solver_1d.hpp"
# include "array_alloc.hpp"

void TEST_tridiag_general();
void TEST_tridiag_ferrari();
void TEST_set_initial_condition_basin();
void TEST_set_initial_condition_periodic();

#endif // UNIT_TESTS_HPP
