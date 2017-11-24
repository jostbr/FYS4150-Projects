
#include "rossby_solver.hpp"

rossby_solver::rossby_solver(double dx, int N, double dt, double T, std::string fileout){
    this->dx = dx;
    this->dt = dt;
    this->T = T;
    this->N = N;
    this->fileout = fileout;

    alloc_array_1D(this->psi_0, N);
    alloc_array_1D(this->zeta_0, N);

    this->outfile.open(fileout);
    this->outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
}


/* Function to retrieve the previously set initial conditions. */
void rossby_solver::get_initial_condition(double* psi_holder, double* zeta_holder) const {
    for (int i = 0; i < this->N; i++){
        psi_holder[i] = this->psi_0[i];
        zeta_holder[i] = this->zeta_0[i];
    }
}


/* Function that writes a time and a psi array to file (on one line). */
void rossby_solver::write_state_to_file(double t, double* psi) {
    this->outfile << std::setw(15) << std::setprecision(8) << t;

    for (int i = 0; i < this->N; i++){
        this->outfile << std::setw(15) << std::setprecision(8) << psi[i];
    }

    this->outfile << std::endl;
}

/* Function for printing the numerical model parameters to screen. */
void rossby_solver::display_model_config() const {
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "dx = " << std::setw(8) << this->dx << std::endl;
    std::cout << "dt = " << std::setw(8) << this->dt << std::endl;
    std::cout << "N = " << std::setw(8) << this->dx << std::endl;
    std::cout << "T = " << std::setw(8) << this->dx << std::endl;
    std::cout << "Output in " << std::setw(8) << this->fileout << std::endl << std::endl;
}


/* Destructor cleaning performing cleanup upon object destruction. */
rossby_solver::~rossby_solver(){
    free_array_1D(this->psi_0);
    free_array_1D(this->zeta_0);
    this->outfile.close();
}
