
#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H

# include <iostream>
# include "ising.h"
# include "mem_alloc.h"

void TEST_get_periodic_index();
void TEST_compute_energy_and_moment();
void TEST_initialize_spin_config_ordered();
void TEST_initialize_spin_config_rng();

#endif // UNIT_TESTS_H
