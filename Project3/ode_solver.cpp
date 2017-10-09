
#include "ode_solver.h"

ODE_solver::ODE_solver(){
    h = 0.0;
    r[0] = 0.0; r[1] = 0.0;
    v[0] = 0.0; v[1] = 0.0;
    a[0] = 0.0; a[1] = 0.0;
}

ODE_solver::~ODE_solver(){

}
