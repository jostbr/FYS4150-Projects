
#include "rossby_solver.hpp"

/* Constructor for a 1D Rossby wave system. Input model parameters for the simulation. */
rossby_solver::rossby_solver(double dx, int N_x, double dt, double T, std::string fileout){
    this->system_dim = "1D";
    this->dx = dx;
    this->N_x = N_x;
    this->dt = dt;
    this->T = T;
    this->fileout = fileout;

    alloc_array_1D(this->psi_0, N_x);
    alloc_array_1D(this->zeta_0, N_x);

    this->outfile.open(fileout);
    this->outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
}


/* Constructor for a 1D Rossby wave system. Input model parameters for the simulation. */
rossby_solver::rossby_solver(double dx, double dy, int N_x, int N_y, double dt, double T, std::string fileout){
    this->system_dim = "2D";
    this->dx = dx;
    this->dy = dy;
    this->N_x = N_x;
    this->N_y = N_y;
    this->dt = dt;
    this->T = T;
    this->fileout = fileout;

    alloc_array_1D(this->psi_0, N_x*N_y);
    alloc_array_1D(this->zeta_0, N_x*N_y);

    this->outfile.open(fileout);
    this->outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
}


/* Function to retrieve the previously set initial conditions. */
void rossby_solver::get_initial_condition(double* psi_holder, double* zeta_holder) const {
    if (this->system_dim.compare("1D") == 0){
        for (int i = 0; i < this->N_x; i++){
            psi_holder[i] = this->psi_0[i];
            zeta_holder[i] = this->zeta_0[i];
        }
    }

    else if (this->system_dim.compare("2D") == 0){
        for (int i = 0; i < this->N_x; i++){
            for (int j = 0; j < this->N_y; j++){
                psi_holder[i*this->N_y + j] = this->psi_0[i*this->N_y + j];
                zeta_holder[i*this->N_y + j] = this->zeta_0[i*this->N_y + j];
            }
        }
    }
}


/* Function that writes a time and a psi array to file (on one line). */
void rossby_solver::write_state_to_file(double t, double* psi){
    this->outfile << std::setw(15) << std::setprecision(8) << t;

    if (this->system_dim.compare("1D") == 0){
        for (int i = 0; i < this->N_x; i++){
            this->outfile << std::setw(15) << std::setprecision(8) << psi[i];
        }

    }

    else if (this->system_dim.compare("2D") == 0){
        for (int j = 0; j < this->N_y; j++){
            for (int i = 0; i < this->N_x; i++){
                this->outfile << std::setw(15) << std::setprecision(8) << psi[i*this->N_y + j];
            }
        }

    }

    this->outfile << std::endl;
}

/* Function for printing the numerical model parameters to screen. */
void rossby_solver::display_model_config() const {
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "dx  = " << std::setw(8) << this->dx << std::endl;

    if (this->system_dim.compare("2D") == 0){
        std::cout << "dy  = " << std::setw(8) << this->dy << std::endl;
    }

    std::cout << "N_x = " << std::setw(8) << this->N_x << std::endl;

    if (this->system_dim.compare("2D") == 0){
        std::cout << "N_y = " << std::setw(8) << this->N_y << std::endl;
    }

    std::cout << "dt  = " << std::setw(8) << this->dt << std::endl;
    std::cout << "T   = " << std::setw(8) << this->T << std::endl;
    std::cout << "\nOutput in " << std::setw(8) << this->fileout << std::endl << std::endl;
}


/* Destructor cleaning performing cleanup upon object destruction. */
rossby_solver::~rossby_solver(){
    free_array_1D(this->psi_0);
    free_array_1D(this->zeta_0);
    this->outfile.close();
}
