
# include <iostream>
# include <cmath>
# include <armadillo>

/* Function that assumes a symmetric matrix A of dimensions NxN. Then finds the
 * maximum value in the matrix excluding values along the diagonal. */
double get_max_non_diag(arma::mat A, int N, int* k, int* l){
    double max_val = 0.0; // Initialize to zero (checking for abs() below)

    /* Only need to loop through upper triangular part of matrix since
     * matrix is assumed symmetric. This also ensures that k < l. */
    for (int i = 0; i < N; i++){
        for (int j = i+1; j < N; j++){      // Ensures i < j --> k < l
            if (fabs(A(i,j)) > fabs(max_val)){
                max_val = fabs(A(i,j));
                *k = i;                     // Store row index for max val
                *l = j;                     // Store column index for max val
            }
        }
    }

    return max_val;
}

/* Function that computes the values for cosine(theta) and sine(theta) needed for
 * the transformation matrix. The values are computed by requiring the transformed
 * element B(k,l) to be zero. This puts a constraint on theta. */
void get_trig_values(arma::mat A, int k, int l, double* cosine, double* sine){
    if (A(k,l) != 0.0){         // More work needed if A(k,l) != 0
        double tau = (A(l,l) - A(k,k))/(2*A(k,l));              // Collecting variables
        double tangent_1 = -tau + sqrt(1 + tau*tau);      // Largest root
        double tangent_2 = -tau - sqrt(1 + tau*tau);      // Smallest root

        *cosine = 1/(sqrt(1 + pow(tangent_2, 2.0)));    // Trig identity for cosine
        *sine = tangent_2*(*cosine);                    // Trig identity for sine
    }

    else {      // Otherwise the values simplify
        *cosine = 1.0;
        *sine = 0.0;
    }
}


/* Function that computes the orthogonal transformation B = S^T*A*S with orthogonal matrix
 * S (a rotation matrix). Due to S being sparse, only some elements in B differ (and thus
 * need to be updated) from those in A. Near the diagonal and around indices k and l. */
void ortho_transform(arma::mat* A, arma::mat* V, int N, int k, int l, double cosine, double sine){
    double a_kk = (*A)(k,k);    // To avoid using overwritten value below
    double a_ll = (*A)(l,l);    // To avoid using overwritten value below
    double a_ik, a_il, v_ik, v_il;

    /* Loop through diagonal of matrix. */
    for (int i = 0; i < N; i++){
        /* Update values next to diagonal. */
        if ((i != k) && (i != l)){
            a_ik = (*A)(i,k);       // To avoid using overwritten value below
            a_il = (*A)(i,l);       // To avoid using overwritten value below
            (*A)(i,k) = a_ik*cosine - a_il*sine;
            (*A)(k,i) = (*A)(i,k);      // Dealing with a symmetric matrix
            (*A)(i,l) = a_il*cosine + a_ik*sine;
            (*A)(l,i) = (*A)(i,l);      // Dealing with a symmetric matrix
        }

        v_ik = (*V)(i,k);
        v_il = (*V)(i,l);
        (*V)(i,k) = v_ik*cosine - v_il*sine;
        (*V)(i,l) = v_il*cosine + v_ik*sine;
    }

    /* Update the two relevant values on the diagonal. */
    (*A)(k,k) = a_kk*cosine*cosine - 2*(*A)(k,l)*cosine*sine
            + a_ll*sine*sine;
    (*A)(l,l) = a_ll*cosine*cosine + 2*(*A)(k,l)*cosine*sine
            + a_kk*sine*sine;
    (*A)(k,l) = 0.0;    // Defined cosine, sine require this to be zero
    (*A)(l,k) = 0.0;    // Defined cosine, sine require this to be zero
}

/* Function that operates as a master function which oversees and executes the
 * jacobi rotation algorithm by calling on the other functions in this file. */
void jacobi_eigen(arma::mat* A, arma::mat* V, int N){
    int k, l;                           // Indices for transf. matrix
    double sine, cosine;                // Trig values for transf. matrix
    double curr_max_nd = 1.0E+10;       // Start with something big

    double eps = 10.0E-10;   // Tolerance for max non-diag element
    int curr_iter = 0;
    int max_iter = N*N*N;

    while ((curr_max_nd > eps) && (curr_iter < max_iter)){
        curr_max_nd = get_max_non_diag(*A, N, &k ,&l);  // Find max non-diag element
        get_trig_values(*A, k, l, &cosine, &sine);      // Find transformation matrix
        ortho_transform(A, V, N, k, l, cosine, sine);   // Execute transformation
        curr_iter++;
    }

    //std::cout << "\nNumber of iterations done: " << curr_iter << std::endl;
}
