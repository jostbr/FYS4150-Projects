
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstdlib>
# include <cmath>
# include <ctime>
# include <string>
# include <armadillo>
# include "assert.h"

/* Functions defined (with description) further below in the program. */
void initialize(double *&, double *&, double *&, double *&, double *&, double *&, int);
void tridiag_general(double*, double*, double*, double*, int, double*);
void tridiag_special(double*, double*, int, double*);
void test_tridiag_general();
void lu_solver();
double source_term(double);
double exact(double);
void time_algorithm(std::string, int, int);
void write_results_to_file(const char*, double*, double*, int);

int main(int argc, char* argv[]){
    for (int i = 1; i < argc; i++){
        int N = atoi(argv[i]);              // Number of points in interior of domain

        std::string fileout = "result";     // Initial part of output filename
        fileout.append(argv[i]);            // Append current problem size to filename
        fileout.append(".txt");

        double *a, *b, *c, *x, *g, *interior;      // See initialize for description of vars.
        double* numerical = new double[N+2];   // To hold numerical solution at all points

        initialize(a, b, c, x, g, interior, N);  // Alloc memory and initialize arrays

        clock_t start_time = clock();
        /* ---------------------------------------------------------------------------------------------- */
        /* Solving the linear system: Choose which algorithm you want by uncommenting one of the two below. */
        //tridiag_general(a, b, c, g, N, interior);         // Get solution in interior with general algorithm
        tridiag_special(b, g, N, interior);            // Get solution in interior with specialized algorithm
        /* ---------------------------------------------------------------------------------------------- */
        clock_t end_time = clock();
        double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;

        numerical[0] = numerical[N+1] = 0.0;                          // Boundary conditions
        for (int i = 1; i < N+1; i++) numerical[i] = interior[i-1];   // Filling full numerical with interior solution

        write_results_to_file(fileout.c_str(), x, numerical, N+2);    // Save results to file

        /* Computing and writing stats to screen. */
        int FLOP = 9*(N - 1) + 1;
        double FLOPS = FLOP/time_used;
        std::cout << "N = " << N << ":" << std::endl;
        std::cout << "--------------------------------" << std::endl;
        std::cout << "Execution time: " << time_used << " seconds" << std::endl;
        std::cout << "FLOPS in algorithm: " << FLOPS << std::endl;
        std::cout << "Writing results to " << fileout << std::endl << std::endl;

        delete[] a; delete[] b; delete[] c; delete[] x; delete[] g;     // Free up memory
        delete[] interior; delete[] numerical;                          // Free up memory
    }

    time_algorithm("general", 3, 5);   // Run timing benchmarks for general algorithm
    time_algorithm("special", 3, 5);   // Run timing benchmarks for special algorithm
    lu_solver();
    test_tridiag_general();     // Unit test on the general algorithm

    return 0;   // End of main function
}


/* Function that allocates memory for diagonals a, b, c, to x-array, right-hand-side g and a an array to hold
 * the solution. The function also initializes a, b, c, x and g with values. */
void initialize(double *& a, double *& b, double *& c, double *& x, double *& g, double *& sol, int N){
    a = new double[N-1];               // To hold lower diagonal
    b = new double[N];                 // To hold main diagonal
    c = new double[N-1];               // To hold upper diagonal
    x = new double[N+2];               // Array to hold all x-values
    g = new double[N];                 // Right side of Av = g
    sol = new double[N];               // To hold num. solution at interior points
    double h = 1.0/((double)(N + 1));  // Step size between each grid point
    double h_squared = h*h;            // To save computing h*h every loop below

    for (int i = 0; i < N-1; i++) a[i] = -1.0;      // Initialize lower diagonal
    for (int i = 0; i < N-1; i++) c[i] = -1.0;      // Initialize upper diagonal
    for (int i = 0; i < N; i++) b[i] = 2.0;         // Initialize main diagonal
    for (int i = 0; i < N+2; i++) x[i] = i*h;       // Compute x-coordinate at all points
    for (int i = 0; i < N; i++) g[i] = h_squared*source_term(x[i+1]);   // Compute right-hand-side
}


/* Function for solving the linear system Ax = y where A is a general tridiagonal matrix (NxN) with
 * lower diag a (length N-1), main diag b (length N) and upper diag c (length N-1). Further y is the
 * known right-hand-side vector and solution is the array to hold the solution x. */
void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution){
    /* STEP 1: Forward substitution. */
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
void tridiag_special(double* b, double* y, int N, double* solution){
    /* STEP 1: Forward substitution. */
    for (int i = 1; i < N; i++){
        b[i] = (i + 2)/((double)(i + 1));       // Eliminating lower diagonal
        y[i] = y[i] + (y[i-1]/b[i-1]);      // Corresponding change to RHS of eq.
    }

    /* STEP 2: Backward substitution. */
    solution[N-1] = y[N-1]/b[N-1];          // Special case for obtaining final element of solution

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] + solution[i+1])/b[i];  // Eliminating upper diag and dividing by main diag
    }
}

/* Function that acts as a unit test for the general tridiag solver. The test is based on using small
 * arbitrary tridiagonal matrix (3x3) with a arbitrary known right-hand-side. We use values for which
 * we've already found a solution by hand. Wwe compare the computed values by tridiag_general() and
 * conclude whether or not the function works as it should. This test function can be helpful during
 * development of the algorithm; to see that any changes made doesn't result in failing the test. */
void test_tridiag_general(){
    int N = 3;
    double* a = new double[N-1];        // Lower diagonal
    double* b = new double[N];          // Main diagonal
    double* c = new double[N-1];        // Upper diagonal
    double* g = new double[N];          // Right hand side of Av = g
    double* solution = new double[N];   // To hold solution from algorithm

    /* Testing with an abitrary matrix for the general algorithm. */
    a[0] = a[1] = c[0] = c[1] = 1.0;
    b[0] = 1.0;
    b[1] = 2.0;
    b[2] = 3.0;
    g[0] = g[1] = g[2] = -1.0;      // Right hand side fortest of general

    std::cout << std::endl << "Running unit test on general tridiagonal algorithm" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    tridiag_general(a, b, c, g, N, solution);   // Call the general algorithm

    /* Examine results and see if results are as expected. Here we chose to use if/else,
     * but we could have used the assert function that auto-terminates program when test fails. */
    if (solution[0] == -1.5 && solution[1] == 0.5 && solution[2] == -0.5){
        std::cout << "Test for general algorithm passed! Behaviour "
                  << "is as expected for the sample matrix." << std::endl << std::endl;
    }

    else {
        std::cout << "Test failed for general algorithm!" << std::endl << std::endl;
    }

    delete[] a; delete[] b; delete[] c; delete[] solution;
}


/* Function that solves the linear algebra problem through LU-decpmosition using the Armadillo
 * library. This functions is designed to compare runtime with the tridiagonal algorithms. */
void lu_solver(){
    std::cout << std::endl << "Benchmarking Armadillo solve algorithm" << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    int upper_exponent = 4;     // 10^upper_exponent is upper limit for N

    for (int k = 1; k < upper_exponent; k++){
        int N = (int)pow(10.0, (double)(k));    // Increase N by factor 10
        int N_full = N + 2;                     // Full domain size
        double h = 1.0/((double)(N_full - 1));  // Step size for current N
        double h_squared = h*h;                 // To ease computaitons for g

        arma::mat A(N, N);      // To hold full matrix
        arma::vec x(N_full);    // To hold x-coordinates
        arma::vec g(N);         // To hold source term on right-hand-side of eq.

        /* Initialize x at all points and g at interior points. */
        for (int i = 0; i < N_full; i++) x[i] = i*h;
        for (int i = 0; i < N; i++) g[i] = h_squared*source_term(x[i+1]);

        /* Initialize matrix with correct elements. */
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                if (i == j){
                    A(i, j) = 2.0;              // Initialize main diagonal
                }

                else if (fabs(i - j) == 1){
                    A(i, j) = -1.0;             // Initialize upper/lower diagonal
                }

                else {
                    A(i, j) = 0.0;      // Initialize zero everywhere else
                }
            }
        }

        arma::vec v;    // Array to hold solution from Armadillo solve algorithm

        double num_runs_per_N = 5;
        double time_total = 0.0;
        clock_t start_time, end_time;

        /* Using Armadillo algorithm many times to find representative average runtime. */
        for (int i = 0; i < num_runs_per_N; i++){
            start_time = clock();           // Start timing
            v = arma::solve(A, g);          // Solve uses LU decomposition internally
            end_time = clock();             // End timing
            double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
            time_total += time_used;

        }

        double time_average = time_total/num_runs_per_N;
        std::cout << "N = " << std::setw(10) << N << "  |  Average execution time: "
                  << std::setw(12) << std::setprecision(8) << time_average << " s" << std::endl;

        /* Testing that Armadillo gives correct solution by
         * writing to file and then plotting with Python. */
        double* numerical = new double[N_full];
        numerical[0] = numerical[N_full-1] = 0.0;   // Boundary conditions

        for (int i = 1; i < N_full-1; i++){
            numerical[i] = v[i-1];              // Fill interior solution to full solution
        }

        write_results_to_file("arma_results.txt", x.memptr(), numerical, N);    // Write solution to file
        delete[] numerical;
    }
}


/* Function returning the value of the source term for a given x. */
double source_term(double x){
    return 100.0*exp(-10.0*x);
}


/* Function returning the ecaxt solution for a given. */
double analytical(double x){
    return 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
}

/* Function that runs the requested algorithm for a series
 * of N-values and then writes the runtime to file. */
void time_algorithm(std::string algorithm, int num_Ns, int num_runs_per_N){
    std::cout << std::endl << "Benchmarking the " << algorithm << " tridiagonal algorithm" << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;

    std::ofstream outfile;                  // File object for output
    std::string time_filename = "timing_";  // Filename for output file
    time_filename.append(algorithm);        // Customize filename
    time_filename.append(".txt");
    outfile.open(time_filename);
    outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    outfile << std::setw(8) << "N" << std::setw(20) << "time" << std::endl;  // Header

    for (int i = 0; i < num_Ns; i++){
        int N = (int)pow(10.0, (double)(i+1));       // Doubling N every iteration
        //double h = 1.0/((double)(N + 1));
        double time_total = 0.0;                    // To keep compute average of num_runs_per_N

        for (int j = 0; j < num_runs_per_N; j++){
            double *a, *b, *c, *x, *g, *solution;       // Need as nput to the algorithms
            initialize(a, b, c, x, g, solution, N);     // Allocate mem. and initialize arrays

            clock_t start_time = clock();           // Start timing

            if (!algorithm.compare("general")){
                tridiag_general(a, b, c, g, N, solution);  // General algorithm
            }

            else if (!algorithm.compare("special")){
                tridiag_special(b, g, N, solution);        // Special algorithm
            }

            clock_t end_time = clock();             // End timing
            double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
            time_total += time_used;    // Increment with current run
            delete[] a; delete[] b; delete[] c; delete[] x; delete[] g; delete[] solution;
        }

        double time_average = time_total/num_runs_per_N;    // Take evareg of num_runs_per_N
        std::cout << "N = " << std::setw(10) << N << "  |  Average execution time: "
                  << std::setw(12) << std::setprecision(8) << time_average << " s" << std::endl;

        /* Write time data to file. */
        outfile << std::setw(10) << N;
        outfile << std::setw(20) << std::setprecision(8) << time_average << std::endl;
    }

    outfile.close();
}


/* Function that writes results x, numerical, analytical and relative error to file filename. */
void write_results_to_file(const char* filename, double* x, double* numerical, int N){
    std::ofstream ofile;    // File object for output file
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "       x:             approx:          exact:       relative error" << std::endl;

    for (int i = 0; i < N; i++) {
        double x_val = x[i];
        double exact_val = analytical(x[i]);
        double rel_error = fabs((exact_val - numerical[i])/exact_val);
        ofile << std::setw(20) << std::setprecision(8) << x_val;
        ofile << std::setw(20) << std::setprecision(8) << numerical[i];
        ofile << std::setw(20) << std::setprecision(8) << exact_val;
        ofile << std::setw(20) << std::setprecision(8) << rel_error << std::endl;
    }

    ofile.close();
}
