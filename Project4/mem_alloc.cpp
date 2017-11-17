
# include "mem_alloc.h"

void alloc_array_2D(int**& array, int num_rows, int num_cols){
    array = new int*[num_rows];

    for(int i = 0; i < num_rows; i++){
        array[i] = new int[num_cols];
    }
}

void free_array_2D(int**& array, int num_rows){
    for(int i = 0; i < num_rows; i++){
        delete [] array[i];
    }

    delete [] array;
}
