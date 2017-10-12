
# include "nbody_solver.h"

nbody_solver::nbody_solver(planet* bodies, int num_bodies){
    this->num_bodies = num_bodies;
    this->bodies = new planet[num_bodies];
    this->ofiles = new std::ofstream[num_bodies];

    for (int i = 0; i < num_bodies; i++){
        this->bodies[i] = bodies[i];
    }
}

void nbody_solver::solve(double h, double t_max, std::string method){
    double frame_write = 100;

    /* Initialize data output file for each bodie. */
    for (int i = 0; i < this->num_bodies; i++){
        std::string filename = this->bodies[i].name;
        filename.append(".txt");
        this->ofiles[i].open(filename.c_str());
        this->ofiles[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        this->ofiles[i] << "t x y" << std::endl;   // Write header to file
        this->write_row_to_file(i, 0.0, this->bodies[i].r[0], this->bodies[i].r[1]);
    }

    if (method.compare("euler") == 0){
        this->euler(h, t_max, frame_write);
    }

    else if (method.compare("verlet") == 0){
        this->verlet(h, t_max, frame_write);
    }

    else {
        std::cout << "\nError: Invalid method! Choose either 'euler' or 'verlet'" << std::endl;
        std::cout << "Terminating program.." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->num_bodies; i++){
        this->ofiles[i].close();   // Close all file objects after time loop
    }
}

void nbody_solver::euler(double h, double t_max, int frame_write){
    planet* bodies_curr = new planet[this->num_bodies];
    int frame = 0;
    double t = 0.0;

    /* Loop until maximum times is reached. */
    while (t <= t_max){
        t += h;     // Increase time by step size
        frame++;    // Count to next frame

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Copy of all planets for current time step
        }

        /* Solve differential equations for all bodies. */
        for (int i = 0; i < this->num_bodies; i++){
            //std::cout << "Time stepping for " << nbodies[i].name << std::endl;

            /* Execute euler time-stepping for both/all directions. */
            for (int dim = 0; dim < cnst::num_dims; dim++){
                bodies_curr[i].a[dim] = bodies_curr[i].compute_total_acc(bodies_curr, this->num_bodies, dim);
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim];
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h*bodies_curr[i].a[dim];
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1]);
                //std::cout << std::endl;
            }
        }
    }

    delete[] bodies_curr;
}

void nbody_solver::verlet(double h, double t_max, int frame_write){
    planet* bodies_curr = new planet[this->num_bodies];
    int frame = 0;
    double h_squared_half = h*h/2.0;
    double h_half = h/2.0;
    double t = 0.0;

    for (int i = 0; i < this->num_bodies; i++){
        bodies_curr[i] = this->bodies[i];        // Copy of all planets for current time step
        bodies_curr[i].a[0] = this->bodies[i].compute_total_acc(this->bodies, this->num_bodies, 0); // Init. x-acc.
        bodies_curr[i].a[1] = this->bodies[i].compute_total_acc(this->bodies, this->num_bodies, 1); // Init. y-acc.
    }

    /* Loop until maximum times is reached. */
    while (t <= t_max){
        t += h;     // Increase time by step size
        frame++;    // Count to next frame

        /* Update positions for all bodies in all directions. */
        for (int i = 0; i < this->num_bodies; i++){
            //std::cout << "Time stepping for " << nbodies[i].name << std::endl;

            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim] +
                        h_squared_half*bodies_curr[i].a[dim];
            }
        }

        /* New loop over bodies and directions to compute the updated values for acceleration
         * and then velocities. Can't do in same loop above since all r_{i+1} must be done for a_{i+1}. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].a[dim] = this->bodies[i].compute_total_acc(
                            this->bodies, this->num_bodies, dim);
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h_half*
                        (this->bodies[i].a[dim] + bodies_curr[i].a[dim]);
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1]);
                //std::cout << std::endl;
            }
        }

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Update current bodies
        }
    }

    delete[] bodies_curr;
}

void nbody_solver::write_row_to_file(int file_index, double t, double x, double y){
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << t;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << x;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << y << std::endl;
}

/* Destructor that deallocates memory from dynamically allocated member variables. */
nbody_solver::~nbody_solver(){
    delete[] this->bodies;
    delete[] this->ofiles;
}

