
# include "nbody_solver.h"

nbody_solver::nbody_solver(){
    num_bodies = 0;
    t_max = 0.0;
    h = 0.0;
}

nbody_solver::nbody_solver(int body_count, double t_end, double step){
    num_bodies = body_count;
    t_max = t_end;
    h = step;
}

void nbody_solver::euler(planet* nbodies){
    std::ofstream ofile[3];
    int num_dims = 2;
    double r_norm_cubed;
    double force_dim_i = 0.0;
    double t = 0.0;
    int frame = 0;

    for (int i = 0; i < num_bodies; i++){
        std::string filename = nbodies[i].name;
        filename.append(".txt");
        ofile[i].open(filename.c_str());
        ofile[i] << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        ofile[i] << "t x y" << std::endl;
        write_row_to_file(ofile[i], t, nbodies[i].r_curr[0], nbodies[i].r_curr[1]);    // Write initial condition to file
    }

    while (t <= t_max){                     // Solve for future
        t += h;
        frame++;

        for (int i = 0; i < num_bodies; i++){   // Solve for both x and y direction
            //std::cout << "Time stepping for " << nbodies[i].name << std::endl;

            for (int dim = 0; dim < num_dims; dim++){      // Solve for all planets
                force_dim_i = nbodies[i].compute_total_force(nbodies, num_bodies, dim);

                nbodies[i].a_curr[dim] = force_dim_i/nbodies[i].mass;             // Compute acceleration in i-direction
                nbodies[i].r_next[dim] = nbodies[i].r_curr[dim] + h*nbodies[i].v_curr[dim];   // Update position in i-direction
                nbodies[i].v_next[dim] = nbodies[i].v_curr[dim] + h*nbodies[i].a_curr[dim];   // Update velocity in i-direction
            }

            if (frame % 100 == 0){
                write_row_to_file(ofile[i], t, nbodies[i].r_next[0], nbodies[i].r_next[1]);    // Write every timestep to file
                std::cout << std::endl;
            }
        }

        /* Set current values (pos, vel) to next values. */
        for (int i = 0; i < num_bodies; i++){
            for (int dim = 0; dim < num_dims; dim++){
                nbodies[i].r_curr[dim] = nbodies[i].r_next[dim];
                nbodies[i].v_curr[dim] = nbodies[i].v_next[dim];
            }
        }

        // Might be something wrong with doing dims loop inside bodies loop concerning the distance computations.
    }

    for (int i = 0; i < num_bodies; i++){
        ofile[i].close();   // Close all file objects
    }
}

void nbody_solver::verlet(){

}

void nbody_solver::write_row_to_file(std::ofstream& ofile, double t, double x, double y){
    ofile << std::setw(20) << std::setprecision(8) << t;
    ofile << std::setw(20) << std::setprecision(8) << x;
    ofile << std::setw(20) << std::setprecision(8) << y << std::endl;
}

nbody_solver::~nbody_solver(){

}

