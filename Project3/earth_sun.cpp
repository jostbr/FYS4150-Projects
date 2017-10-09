
# include <iostream>
# include <iomanip>
# include <fstream>
# include <string>
# include <cmath>

void euler_solver(double* r, double* v, double h, double t_max, std::ofstream& ofile);
void verlet_solver(double* r_0, double* v_0, double h, double t_max, std::ofstream& ofile);
void write_row_to_file(std::ofstream& ofile, double t, double x, double y);

void two_body(){
    std::string filename = "results.txt";
    std::ofstream ofile;    // File object for output file
    ofile.open(filename);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "t x y" << std::endl;

    double pi = acos(-1.0);
    double t_max = 2.0;               // Upper time in years
    double h = 0.0001;         // Step size in years

    double r_0[2];            // Position x, y in AU
    double v_0[2];            // Velcoity v_x, v_y in AU/yr

    double t = 0.0;
    r_0[0] = 1.0; r_0[1] = 0.0;
    v_0[0] = 0.0; v_0[1] = 2.0*pi;

    //euler_solver(r_0, v_0, h, t_max, ofile);
    verlet_solver(r_0, v_0, h, t_max, ofile);

    ofile.close();
}

/* Function Newtons 2nd law with gravity as the force using the euler method. */
void euler_solver(double* r_0, double* v_0, double h, double t_max, std::ofstream& ofile){
    int num_dims = 2;
    double four_pi_sq = 4*acos(-1.0)*acos(-1.0);
    double r_norm_cubed;
    double t = 0.0;
    double r[2], v[2], a[2];

    write_row_to_file(ofile, t, r_0[0], r_0[1]);    // Write initial condition to file

    r[0] = r_0[0]; r[1] = r_0[1];
    v[0] = v_0[0]; v[1] = v_0[1];

    while (t <= t_max){
        r_norm_cubed = pow(sqrt(r[0]*r[0] + r[1]*r[1]), 3.0);

        for (int i = 0; i < num_dims; i++){
            a[i] = -four_pi_sq*r[i]/r_norm_cubed;   // Compute acceleration in i-direction
            r[i] = r[i] + h*v[i];                   // Update position in i-direction
            v[i] = v[i] + h*a[i];                   // Update velocity in i-direction
        }

        t += h;

        write_row_to_file(ofile, t, r[0], r[1]);    // Write every timestep to foile
    }
}


/* Function Newtons 2nd law with gravity as the force using the velocity verlet method. */
void verlet_solver(double* r_0, double* v_0, double h, double t_max, std::ofstream& ofile){
    int num_dims = 2;
    double four_pi_sq = 4*acos(-1.0)*acos(-1.0);
    double h_squared = h*h;
    double h_half = h/2;
    double r_norm_cubed;
    double t = 0.0;
    double r[2], v[2], a_prev[2], a_next[2];

    write_row_to_file(ofile, t, r_0[0], r_0[1]);    // Write initial condition to file

    r[0] = r_0[0]; r[1] = r_0[1];
    v[0] = v_0[0]; v[1] = v_0[1];

    while (t <= t_max){
        r_norm_cubed = pow(sqrt(r[0]*r[0] + r[1]*r[1]), 3.0);

        for (int dim = 0; dim < num_dims; dim++){
            a_prev[dim] = -four_pi_sq*r[dim]/r_norm_cubed;   // Compute acc. in dim-direction for curr time step
            r[dim] = r[dim] + h*v[dim] + (h_squared/2.0)*a_prev[dim];   // Update position in dim-direction

            a_next[dim] = -four_pi_sq*r[dim]/r_norm_cubed;   // Compute acc. in dim-direction for next time step

            v[dim] = v[dim] + h_half*(a_next[dim] + a_prev[dim]);   // Update velocity in dim-direction
        }

        t += h;

        write_row_to_file(ofile, t, r[0], r[1]);    // Write every timestep to foile
    }
}

void write_row_to_file(std::ofstream& ofile, double t, double x, double y){
    ofile << std::setw(20) << std::setprecision(8) << t;
    ofile << std::setw(20) << std::setprecision(8) << x;
    ofile << std::setw(20) << std::setprecision(8) << y << std::endl;
}
