
# include <iostream>
# include <cstdlib>
# include <cmath>

typedef double (*func)(double arg);  // Type for func to differentiate

/* List of functions located further below. */
double arctan(double x);
double dfdx_2c(func f, double x, double h);
double dfdx_3c(func f, double x, double h);

/* Main function doing stuff. */
int main(int argc, char* argv[]){
    double h_start = strtod(argv[1], NULL);
    double h_final = strtod(argv[2], NULL);
    int num_h_values = atoi(argv[3]);

    double h_step = (h_start - h_final)/(num_h_values - 1);
    double x = sqrt(2);     // Where to compute the derivative
    double* h_values = new double[4];

    for (int i = 0; i < num_h_values; i++){
        h_values[i] = h_start - i*h_step;
        std::cout << h_values[i] << std::endl;
    }

    /*double val_1 = dfdx_2c(arctan, x, h);
    double val_2 = dfdx_3c(arctan, x, h);
    std::cout << val_1 << std::endl;
    std::cout << val_2 << std::endl;*/
    return 0;
}

/* Function returning the inverse tangent. */
double arctan(double x){
    return atan(x);
}

/* Function computing first order derivative
using a finite difference forward 2-point formula. */
double dfdx_2c(func f, double x, double h){
    return (f(x + h) - f(x))/h;
}

/* Function computing first order derivative
using a finite difference centered 3-point formula. */
double dfdx_3c(func f, double x, double h){
    return (f(x + h) - f(x - h))/(2*h);
}
