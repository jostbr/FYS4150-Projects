
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstdlib>
# include <cmath>
# include <ctime>
# include <sstream>

/* Functions defined further below in the program. */
void tridiag_general(double*, double*, double*, double*, int, double*);
void tridiag_specialized(double*, int, double*);
double source_term(double);
double exact(double);
void write_results_to_file(const char*, double*, double*, double*, int);

int main(int argc, char* argv[]){
	int upper_exponent = atoi(argv[1]);
    int N;

    for (int i = 1; i <= upper_exponent; i++){
        N = pow(10, i);

        std::string fileout = "result";
        std::ostringstream convert;         // Stream used for the conversion
        convert << N;                       // Insert the text rep of N in the stream
        std::string result = convert.str(); // Set to the contents of the stream
        fileout.append(result);
        fileout.append(".txt");
        std::cout << fileout << std::endl;
    
    	double h = 1.0/((double)(N + 1));
    	double* x = new double[N];

    	double* a = new double[N-1];         // Lower diagonal
    	double* b = new double[N];         // Main diagonal
    	double* c = new double[N-1];         // Upper diagonal
    	double* g = new double[N];         // Right side of diff. eq.
    	double* interior = new double[N];  // To hold solutionution
    	double* numerical = new double[N+2];   // To hold solutionution
    	double* analytical = new double[N+2];  // To hold analytical solusolution
        double h_squared = h*h;              // To save computing h*h every loop

        for (int i = 0; i < N-1; i++) a[i] = -1.0;
    	for (int i = 0; i < N-1; i++) c[i] = -1.0;
        for (int i = 0; i < N; i++) b[i] = 2.0;

        for (int i = 0; i < N+2; i++){
    		x[i] = i*h;
    		analytical[i] = exact(x[i]);
    	}

        for (int i = 0; i < N; i++){
            g[i] = h_squared*source_term(x[i+1]);
        }

    	numerical[0] = numerical[N+1] = 0.0;              // Boundary conditions
        clock_t start_time = clock();
    	tridiag_general(a, b, c, g, N, interior);     // Solution for interior
        //tridiag_specialized(g, N, interior);
        clock_t end_time = clock();
        double time_used = (double)(end_time - start_time)/CLOCKS_PER_SEC;
        std::cout << std::setw(15) << std::setprecision(10) << "Time :" << time_used << "s" << std::endl;

    	for (int i = 1; i < N+1; i++){
    		numerical[i] = interior[i-1];
    	}

        write_results_to_file(fileout.c_str(), x, numerical, analytical, N+2);
    }

	return 0;
}

void tridiag_general(double* a, double* b, double* c, double* y, int N, double* solution){
	for (int i = 1; i < N; i++){
		b[i] = b[i] - a[i-1]*a[i-1]/b[i-1];
		y[i] = y[i] - (a[i-1]/b[i-1])*y[i-1];
	}

	solution[N-1] = y[N-1]/b[N-1];

	for (int i = N-2; i >= 0; i--){
		solution[i] = (y[i] - c[i]*solution[i+1])/b[i];
	}
}

void tridiag_specialized(double* y, int N, double* solution){
    double* b = new double[N];
    b[0] = 2;

    for (int i = 1; i < N; i++){
        b[i] = (i + 2)/((double)(i + 1));
        y[i] = y[i] + (y[i-1]/b[i-1]);
    }

    solution[N-1] = y[N-1]/b[N-1];

    for (int i = N-2; i >= 0; i--){
        solution[i] = (y[i] + solution[i+1])/b[i];
    }
}


double source_term(double x){
    return 100.0*exp(-10.0*x);
}


double exact(double x){
    return 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
}


void write_results_to_file(const char* filename, double* x, double* numerical, double* analytical, int N){
    /* Write results to output file. */
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "       x:             approx:          exact:       relative error" << std::endl;

    for (int i = 0; i < N; i++) {
        double x_val = x[i];
        double rel_error = fabs((exact(x_val)-numerical[i])/exact(x_val));
        ofile << std::setw(20) << std::setprecision(8) << x_val;
        ofile << std::setw(20) << std::setprecision(8) << numerical[i];
        ofile << std::setw(20) << std::setprecision(8) << exact(x_val);
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