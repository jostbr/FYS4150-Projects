
#ifndef ARRAY_ALLOC_HPP
#define ARRAY_ALLOC_HPP

void alloc_array_1D(double*& array, int num_elements);
void free_array_1D(double*& array);
void alloc_array_2D(double**& array, int num_rows, int num_cols);
void free_array_2D(double**& array, int num_rows);

#endif // ARRAY_ALLOC_HPP
