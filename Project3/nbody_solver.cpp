
# include "nbody_solver.h"

nbody_solver::nbody_solver(planet* bodies, int num_bodies){
    this->num_bodies = num_bodies;
    this->bodies = new planet[num_bodies];

    for (int i = 0; i < num_bodies; i++){
        this->bodies[i] = bodies[i];
    }
}

void nbody_solver::euler(double h, double t_max, std::string method){
    planet* bodies_curr = new planet[this->num_bodies];
    std::ofstream* ofiles = new std::ofstream[this->num_bodies];
    int num_dims = 2;
    double t = 0.0;
    double a_curr[2];
    int frame = 0;

    this->initialize_output_files(&ofiles);

    /* Loop until maximum times is reached. */
    while (t <= t_max){
        t += h;     // Increase time by step size
        frame++;    // Count to next frame

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Make copy of all planets for current time step
        }

        for (int i = 0; i < this->num_bodies; i++){   // Solve for both x and y direction
            //std::cout << "Time stepping for " << nbodies[i].name << std::endl;
            for (int dim = 0; dim < num_dims; dim++){      // Solve for all planets
                a_curr[dim] = bodies_curr[i].compute_acceleration(bodies_curr, this->num_bodies, dim); // Compute acceleration in dim-direction
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim];   // Update position in i-direction
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h*a_curr[dim];   // Update velocity in i-direction
            }

            if (frame % 100 == 0){
                this->write_row_to_file(ofiles[i], t, this->bodies[i].r[0], this->bodies[i].r[1]);    // Write every timestep to file
                //std::cout << std::endl;
            }
        }
    }

    for (int i = 0; i < this->num_bodies; i++){
        ofiles[i].close();   // Close all file objects after time loop
    }

    delete[] bodies_curr; delete[] ofiles;
}

void nbody_solver::verlet(){

}

void nbody_solver::initialize_output_files(std::ofstream** ofiles){
    for (int i = 0; i < this->num_bodies; i++){
        std::string filename = this->bodies[i].name;
        filename.append(".txt");
        (*ofiles)[i].open(filename.c_str());
        (*ofiles)[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        (*ofiles)[i] << "t x y" << std::endl;   // Write header to file
        this->write_row_to_file((*ofiles)[i], 0.0, this->bodies[i].r[0], this->bodies[i].r[1]);
    }
}

void nbody_solver::write_row_to_file(std::ofstream& ofile, double t, double x, double y){
    ofile << std::setw(20) << std::setprecision(8) << t;
    ofile << std::setw(20) << std::setprecision(8) << x;
    ofile << std::setw(20) << std::setprecision(8) << y << std::endl;
}

nbody_solver::~nbody_solver(){

}

