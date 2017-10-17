
# include "nbody_solver.h"

/* Overload constructor (no ordinary exists) function with arguments. */
nbody_solver::nbody_solver(planet* bodies, int n, bool implicit_sun){
    this->num_bodies = n;
    this->fixed_sun = implicit_sun;
    this->bodies = new planet[num_bodies];
    this->ofiles = new std::ofstream[num_bodies];

    for (int i = 0; i < num_bodies; i++){
        this->bodies[i] = bodies[i];
    }
}

/* Manager function that supervises the n-body solution process. This includes ccalling a solution
 * algorithm, making sure parameter arguments are reasonable and handling ouput data for each body. */
void nbody_solver::solve(double h, double t_max, double t_write, std::string method){
    int frame_write = (int)(t_write/h + 0.5);        // Write to file interval
    int total_frames = (int)(t_max/t_write + 0.5);

    std::cout << "\nSolving n-body problem with following parameters:" << std::endl;
    std::cout << "=================================================" << std::endl;
    std::cout << "Using algorithm :           " << method << std::endl;
    std::cout << "Number of bodies:           " << this->num_bodies << std::endl;
    std::cout << "Simulation time [years]:    " << t_max << std::endl;
    std::cout << "Time step, h [years]:       " << std::setprecision(2) << h << std::endl;
    std::cout << "Total time steps:           " << t_max/h << std::endl;
    std::cout << "Output every [time step]:   " << frame_write << std::endl;
    std::cout << "Total output [time steps]:  " << total_frames << std::endl << std::endl;

    /* If arguments passed would result in a "dangerously" large output files. */
    if (total_frames > 50000){
        std::cout << "\nError: Invalid choice of t_max and t_write!" << std::endl;
        std::cout << "Lines written to file can't exceed t_max/t_write = " << total_frames << std::endl;
        std::cout << "Terminating program...\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Initialize data output file for each body. */
    for (int i = 0; i < this->num_bodies; i++){
        std::string filename = this->bodies[i].name;
        filename.append(".txt");
        this->ofiles[i].open(filename.c_str());
        this->ofiles[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        this->ofiles[i] << "t x y z" << std::endl;   // Write header to file
        this->write_row_to_file(i, 0.0, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
    }

    clock_t t_0 = clock();

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

    clock_t t_1 = clock();

    double time_used = (double)(t_1 - t_0)/CLOCKS_PER_SEC;
    std::cout << "\nTime used by " << method << " method: " << std::setprecision(8)
              << time_used << " seconds\n" << std::endl;

    for (int i = 0; i < this->num_bodies; i++){
        this->ofiles[i].close();   // Close all file objects after time loop
    }
}

/* Function that solves the n-body problem using the "Forward Euler" algorithm for time-stepping
 * the position and velocity in the x- and y-direction (can be extended to include z as well). The
 * only assumed knowlede about the body objects used in this function, is that they have a member
 * function compute_total_acc() that computes and returns the total acceleration on a body. */
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
                bodies_curr[i].a[dim] = this->compute_total_acc(bodies_curr[i], bodies_curr, dim);
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim];
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h*bodies_curr[i].a[dim];
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
                //std::cout << std::endl;
            }
        }
    }

    delete[] bodies_curr;
}

/* Function that solves the n-body problem using the "Velocity Verlet" algorithm for time-stepping
 * the position and velocity in the x- and y-direction (can be extended to include z as well). The
 * only assumed knowlede about the body objects used in this function, is that they have a member
 * function compute_total_acc() that computes and returns the total acceleration on a body. */
void nbody_solver::verlet(double h, double t_max, int frame_write){
    planet* bodies_curr = new planet[this->num_bodies];
    int frame = 0;                      // Counting time steps for writing to file
    double h_squared_half = h*h/2.0;
    double h_half = h/2.0;
    double t = 0.0;

    /* Make copy of all bodies to track current and next time step. Also compute init. acceleration. */
    for (int i = 0; i < this->num_bodies; i++){
        bodies_curr[i] = this->bodies[i];

        for (int dim = 0; dim < cnst::num_dims; dim++){
            bodies_curr[i].a[dim] = this->compute_total_acc(this->bodies[i], this->bodies, dim);
        }
    }

    /* Loop until maximum times is reached. */
    while (t <= t_max){
        t += h;     // Increase time by step size
        frame++;    // Count to next frame

        /* Update positions for all bodies in all directions. */
        for (int i = 0; i < this->num_bodies; i++){
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].r[dim] = bodies_curr[i].r[dim] + h*bodies_curr[i].v[dim] +
                        h_squared_half*bodies_curr[i].a[dim];
            }
        }

        /* New loop over bodies and directions to compute the updated values for acceleration
         * and then velocities. Can't do in same loop above since all r_{i+1} must be done for a_{i+1}. */
        for (int i = 0; i < this->num_bodies; i++){
            //std::cout << "Time stepping for " << this->bodies[i].name << std::endl;
            for (int dim = 0; dim < cnst::num_dims; dim++){
                this->bodies[i].a[dim] = this->compute_total_acc(this->bodies[i], this->bodies, dim);
                this->bodies[i].v[dim] = bodies_curr[i].v[dim] + h_half*
                        (this->bodies[i].a[dim] + bodies_curr[i].a[dim]);
            }

            if (frame % frame_write == 0){
                this->write_row_to_file(i, t, this->bodies[i].r[0], this->bodies[i].r[1], this->bodies[i].r[2]);
            }

            //std::cout << std::endl;
        }

        for (int i = 0; i < this->num_bodies; i++){
            bodies_curr[i] = this->bodies[i];        // Update current bodies
        }
    }

    delete[] bodies_curr;
}

/* Function that, through calls to planet::compute_accel(), computes total acceleration of subject
 * (in dim-direction) due to the sum of gravitational forces from all planets in array objects. */
double nbody_solver::compute_total_acc(planet subject, planet* objects, int dim) const {
    double total_accel = 0.0;

    /* Simulation is being run with a fixed sun at the origin. */
    if (this->fixed_sun == true){
        std::cout << "Using fixed sun!" << std::endl;
        planet sun("sun", cnst::mass_sun, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        total_accel += subject.compute_acc(sun, dim);     // Include force from sun separately
    }

    for (int i = 0; i < this->num_bodies; i++){
        if (objects[i].name.compare(subject.name) != 0){    // No force on itself
            //if (dim == 0) std::cout << "While using force from " << planets[i].name << std::endl;
            total_accel += subject.compute_acc(objects[i], dim);
        }
    }

    return total_accel;  // Divide by this->mass to get acceleration
}

/* Function that writes a row of values (t, x, y) to data member output file
 * ofiles[file_index] where file_index is unique per body i the simulation. */
void nbody_solver::write_row_to_file(int file_index, double t, double x, double y, double z){
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << t;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << x;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << y;
    this->ofiles[file_index] << std::setw(20) << std::setprecision(8) << z << std::endl;
}

/* Destructor that deallocates memory from dynamically allocated member variables. */
nbody_solver::~nbody_solver(){
    delete[] this->bodies;
    delete[] this->ofiles;
}

