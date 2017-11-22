
# include "array_alloc.hpp"

/* Function that dynamically allocates a 1D array with specified number of elements. */
void alloc_array_1D(double*& array, int num_elements){
    array = new double[num_elements];
}

/* Function that deallocates a dynamically allocated 1D array. */
void free_array_1D(double*& array){
    delete [] array;
}

/* Function that dynamically allocates a 2D array with specified number of rows and columns. */
void alloc_array_2D(double**& array, int num_rows, int num_cols){
    array = new double*[num_rows];

    for (int i = 0; i < num_rows; i++){
        array[i] = new double[num_cols];
    }
}

/* Function that deallocates a dynamically allocated 2D array. */
void free_array_2D(double**& array, int num_rows){
    for (int i = 0; i < num_rows; i++){
        delete [] array[i];
    }

    delete [] array;
}
