#ifndef RIEMANN_H_
#define RIEMANN_H_

#include "fv.h"

typedef struct riemann_sol
{
    double velocity_star;
    double pressure_star;
    double rho_s_left;
    double rho_s_right;

    int left_state; //0 = shock
    int right_state; //0 = shock

    double SHL;
    double STL;
    double SHR;
    double STR;
} riemann_sol;

riemann_sol solve(cell_state left, cell_state right);

cell_state solve_new_state(riemann_sol sol,  cell_state left, cell_state right);

#endif