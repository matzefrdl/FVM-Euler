#ifndef DISC_H_
#define DISC_H_

#include "riemann_solver.h"

typedef struct cell_grid
{
    cell_state_2d** cells;
    int width;
    int height;
    double dx;
} cell_grid;

void run_timestep(cell_grid* grid, double dt);
cell_grid* create_grid(int width, int height, double dx);
void free_grid(cell_grid* grid);

cell_state_2d compute(cell_state_2d current,
                         cell_state_2d top,
                         cell_state_2d right,
                         cell_state_2d bottom,
                         cell_state_2d left, double dx, double dt);

#endif