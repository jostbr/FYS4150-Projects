
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstdlib>
# include <cmath>

std::ofstream ofile;

void tridiag_gen(double*, double*, double*, double*, int, double**);
double exact(double);

int main(int argc, char* argv[]){
    const char* fileout = "output.txt";
	int N = atoi(argv[1]);
	double h = 1.0/((double)(N - 1));
	double* x = new double[N];

	double* a = new double[N-3];	// Lower diagonal
	double* b = new double[N-2];		// Main diagonal
	double* c = new double[N-3];	// Upper diagonal
	double* y = new double[N-2];		// Right side of diff. eq.
	double* sol = new double[N-2];		// To hold solution
	double* numerical = new double[N];		// To hold solution
	double* analytical = new double[N];

	for (int i = 0; i < N-2; i++){
		if (i < N-3){
			a[i] = -1.0;
			c[i] = -1.0;
		}

		b[i] = 2.0;
        y[i] = h*h*100.0*exp(-10*x[i]);
    }

    for (int i = 0; i < N; i++){
		x[i] = i*h;
		analytical[i] = exact(x[i]);
	}

	numerical[0] = numerical[N] = 0.0;

	tridiag_gen(a, b, c, y, N-2, &sol);

	int j = 0;

	for (int i = 1; i <= N-1; i++){
		numerical[i] = sol[j];
		j += 1;
	}

	ofile.open(fileout);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    //      ofile << "       x:             approx:          exact:       relative error" << endl;
    for (int i = 1; i < N-1; i++) {
        double xval = x[i];
        double RelativeError = fabs((exact(xval)-numerical[i])/exact(xval));
        ofile << std::setw(15) << std::setprecision(8) << xval;
        ofile << std::setw(15) << std::setprecision(8) << numerical[i];
        ofile << std::setw(15) << std::setprecision(8) << exact(xval);
        ofile << std::setw(15) << std::setprecision(8) << log10(RelativeError) << std::endl;
    }

    ofile.close();

	return 0;
}

void tridiag_gen(double* a, double* b, double* c, double* y, int N, double** sol){
	double* b_tilde = new double[N-1];
	double* y_tilde = new double[N-1];

	b_tilde[0] = b[1] - a[0]*a[0]/b[0];
	y_tilde[0] = y[1] - (a[0]/b[0])*y[0];

	for (int i = 1; i < N-1; i++){
		b_tilde[i] = b[i] - a[i-1]*c[i-1]/b_tilde[i-1];
		y_tilde[i] = y[i] - (a[i-1]/b_tilde[i-1])*y_tilde[i-1];
	}

	(*sol)[N-1] = y_tilde[N-2]/b_tilde[N-2];

	for (int i = N-2; i > 0; i--){
		(*sol)[i] = (y[i] - c[i]*(*sol)[i+1])/b_tilde[i];
	}

	(*sol)[0] = (y[0] - c[0]*(*sol)[2])/b[0];
}

double exact(double x){
    return 1 - (1 - exp(-10.0))*x - exp(-10.0*x);
}
