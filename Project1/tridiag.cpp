
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstdlib>
# include <cmath>
# include <ctime>

/* Functions defined (with description) further below in the program. */
void tridiag_general(double*, double*, double*, double*, int, double*);
void tridiag_specialized(double*, int, double*);
double source_term(double);
double exact(double);
void write_results_to_file(const char*, double*, double*, double*, int);

int main(int argc, char* argv[]){
    for (int i = 1; i < argc; i++){
        int N = atoi(argv[i]);              // Number of points in interior of domain

        std::string fileout = "result";     // Initial part of output filename
        fileout.append(argv[i]);            // Append current problem size to filename
        fileout.append(".txt");

        double h = 1.0/((double)(N + 1));   // Step size between each grid point
        double* x = new double[N];          // Array to hold all x-values

        double* a = new double[N-1];           // To hold lower diagonal
        double* b = new double[N];             // To hold main diagonal
        double* c = new double[N-1];           // To hold upper diagonal
        double* g = new double[N];             // Right side of Av = g
        double* interior = new double[N];      // To hold num. solution at interior points
        double* numerical = new double[N+2];   // To hold num. solution at all points
        double* analytical = new double[N+2];  // To hold analytical solution at all points
        double h_squared = h*h;                // To save computing h*h every loop below

        for (int i = 0; i < N-1; i++) a[i] = -1.0;      // Initialize lower diagonal
        for (int i = 0; i < N-1; i++) c[i] = -1.0;      // Initialize upper diagonal
        for (int i = 0; i < N; i++) b[i] = 2.0;         // Initialize main diagonal

        for (int i = 0; i < N+2; i++){
            x[i] = i*h;                     // Compute x-coordinate at all points
            analytical[i] = exact(x[i]);    // Compute analytical solution at all points
        }

        for (int i = 0; i < N; i++){
            g[i] = h_squared*source_term(x[i+1]);   // Compute right-hand-side of Av = g
        }

        numerical[0] = numerical[N+1] = 0.0;              // Boundary conditions
        clock_t start_time = clock();

        /* ---------------------------------------------------------------------------------------------- */
        /* Solving the linear system: Choose which algorithm you want by uncommenting one of the two below. */
        //tridiag_general(a, b, c, g, N, interior);         // Get solution in interior with general algorithm
        tridiag_specialized(g, N, interior);            // Get solution in interior with specialized algorithm
        /* ---------------------------------------------------------------------------------------------- */

        clock_t end_time = clock();
        double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;

        for (int i = 1; i < N+1; i++){
            numerical[i] = interior[i-1];   // Filling full numerical with interior solution
        }

        write_results_to_file(fileout.c_str(), x, numerical, analytical, N+2);

        /* Computing and writing stats to screen. */
        int FLO = 9*(N - 1) + 1;
        double FLOPS = FLO/time_used;
        std::cout << "N = " << N << ":" << std::endl;
        std::cout << "--------------------------------" << std::endl;
        std::cout << "Execution time: " << time_used << " seconds" << std::endl;
        std::cout << "FLOPS in algorithm: " << FLOPS << std::endl;
        std::cout << "Writing results to " << fileout << std::endl << std::endl;
    }

    return 0;   // End of main function
}


/* Function for solving the linear system Ax = y where A is a general tridiagonal matrix (NxN) with
 * lower diag a (length N-1), main diag b (length N) and upper diag c (length N-1). Further y is the
 * known right-hand-side vector and solution is the array to hold the solution x. */
void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution){
    /* STEP 1: Forward substitution and decomposition. */
    for (int i = 1; i < N; i++){
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1];     // Eliminating lower diagonal
        y[i] = y[i] - (a[i-1]/b[i-1])*y[i-1];   // Corresponding change to RHS of eq.
    }

    /* STEP 2: Backward substitution. */
    solution[N-1] = y[N-1]/b[N-1];  // Special case for obtaining final element of solution

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] - c[i]*solution[i+1])/b[i];  // Eliminating upper diag and dividing by main diag
    }
}

/* Function for solving the linear system Ax = y where A is a special tridiagonal matrix (NxN) with
 * all elements along lower and upper diag equal to -1, while the main diag has all values equal to
 * 2. Further y is the known right-hand-side vector and solution is the array to hold the solution x. */
void tridiag_specialized(double* y, int N, double* solution){
    double* b = new double[N];  // To hold substituted values along main diag
    b[0] = 2;                   // First substituted lement is "itself" (2)

    for (int i = 1; i < N; i++){
        b[i] = (i + 2)/((double)(i + 1));
        y[i] = y[i] + (y[i-1]/b[i-1]);
    }

    solution[N-1] = y[N-1]/b[N-1];

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] + solution[i+1])/b[i];
    }
}

/* Function returning the value of the source term for a given x. */
double source_term(double x){
    return 100.0*exp(-10.0*x);
}

/* Function returning the ecaxt solution for a given. */
double exact(double x){
    return 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
}

/* Function that writes results x, numerical, analytical and reslative error to file filename. */
void write_results_to_file(const char* filename, double* x, double* numerical, double* analytical, int N){
    /* Write results to output file. */
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "       x:             approx:          exact:       relative error" << std::endl;

    for (int i = 0; i < N; i++) {
        double x_val = x[i];
        double rel_error = fabs((analytical[i]-numerical[i])/analytical[i]);
        ofile << std::setw(20) << std::setprecision(8) << x_val;
        ofile << std::setw(20) << std::setprecision(8) << numerical[i];
        ofile << std::setw(20) << std::setprecision(8) << analytical[i];
        ofile << std::setw(20) << std::setprecision(8) << log10(rel_error) << std::endl;
    }

    ofile.close();
}

/*void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution){
    double* b = new double[N];
    double* y = new double[N];

    b[0] = b[0];
    y[0] = y[0];

    for (int i = 1; i < N; i++){
        b[i] = b[i] - a[i-1]*a[i-1]/b[i-1];
        y[i] = y[i] - (a[i-1]/b[i-1])*y[i-1];
    }

    solution[N-1] = y[N-1]/b[N-1];

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] - c[i]*solution[i+1])/b[i];
    }
}*/
