#include "discrete.h"
#include "math.h"
#include <stdlib.h>

cell_state flatten_x(cell_state_2d state2d)
{
    cell_state state = {0};
    state.pressure = state2d.pressure;
    state.density = state2d.density;
    state.velocity = state2d.velocity.x;
    return state;
}
cell_state flatten_y(cell_state_2d state2d)
{
    cell_state state = {0};
    state.pressure = state2d.pressure;
    state.density = state2d.density;
    state.velocity = state2d.velocity.y;
    return state;
}

typedef struct time_derivative
{
    double drhodt;
    vec2 dmdt;
    double dEdt;
} time_derivative;

time_derivative compute_dt(cell_state left, cell_state right, double transverse_l, double transverse_r)
{
    double gamma = 1.4;
    riemann_sol sol = solve(left, right);

    cell_state at_interface = solve_new_state(sol, left, right);

    //fake x notation
    double drhodt = -at_interface.density * at_interface.velocity;
    double drhoudt = -at_interface.density * at_interface.velocity * at_interface.velocity - at_interface.pressure;
    
    double transverse = transverse_r;

    if(at_interface.velocity > 0)
    {
        transverse = transverse_l;
    }

    double drhovdt = -at_interface.density * at_interface.velocity * transverse; 

    vec2 dmdt = {drhoudt, drhovdt};

    time_derivative deriv = {0};
    deriv.drhodt = drhodt;
    deriv.dmdt = dmdt;
    double u2 = at_interface.velocity * at_interface.velocity;
    double v2 = transverse*transverse;
    double E = 0.5 * at_interface.density * (u2+v2) + at_interface.pressure / (gamma-1);
    deriv.dEdt = -at_interface.velocity * (E + at_interface.pressure);
    return deriv;
}

cell_state_2d compute(cell_state_2d current,
                         cell_state_2d top,
                         cell_state_2d right,
                         cell_state_2d bottom,
                         cell_state_2d left,
                         double dx, double dt)
{
    double gamma = 1.4;

    cell_state_2d new = {0};
    new.pressure = current.pressure;
    new.velocity = current.velocity;
    new.momenta = current.momenta;
    new.energy = current.energy;
    new.density = current.density;

    //top (y - step)
    cell_state currenty = flatten_y(current);
    cell_state topy = flatten_y(top);

    time_derivative deriv = compute_dt(currenty, topy, current.velocity.x, top.velocity.x);
    new.energy += deriv.dEdt * dt / dx;
    new.density += deriv.drhodt * dt / dx;

    new.momenta.y += deriv.dmdt.x * dt / dx;
    new.momenta.x += deriv.dmdt.y * dt / dx;
    //bottom 
    cell_state bottomy = flatten_y(bottom);

    deriv = compute_dt(bottomy, currenty, bottom.velocity.x, current.velocity.x);
    new.energy -= deriv.dEdt * dt / dx;
    new.density -= deriv.drhodt * dt / dx;
    
    new.momenta.y -= deriv.dmdt.x * dt / dx;
    new.momenta.x -= deriv.dmdt.y * dt / dx;
    //subtract because of \int f' = F(b) - F(a)

    //fix old variables
    new.velocity.x = new.momenta.x / new.density;
    new.velocity.y = new.momenta.y / new.density;
    new.pressure = (gamma - 1) * (new.energy - 0.5 * new.density * (new.velocity.x*new.velocity.x+new.velocity.y*new.velocity.y));

    current = new;

    //x step
    cell_state currentx = flatten_x(current);
    cell_state rightx = flatten_x(right);

    deriv = compute_dt(currentx, rightx, current.velocity.y, right.velocity.y);
    new.energy += deriv.dEdt * dt / dx;
    new.density += deriv.drhodt * dt / dx;

    new.momenta.x += deriv.dmdt.x * dt / dx;
    new.momenta.y += deriv.dmdt.y * dt / dx;
    
    //left 
    cell_state leftx = flatten_x(left);

    deriv = compute_dt(leftx, currentx, left.velocity.y, current.velocity.y);
    new.energy -= deriv.dEdt * dt / dx;
    new.density -= deriv.drhodt * dt / dx;
    
    new.momenta.x -= deriv.dmdt.x * dt / dx;
    new.momenta.y -= deriv.dmdt.y * dt / dx;
    //subtract because of \int f' = F(b) - F(a)

    //fix old variables
    new.velocity.x = new.momenta.x / new.density;
    new.velocity.y = new.momenta.y / new.density;
    new.pressure = (gamma - 1) * (new.energy - 0.5 * new.density * (new.velocity.x*new.velocity.x+new.velocity.y*new.velocity.y));

    current = new;

    return current;
}

int is_wall(int x, int y, cell_grid* grid)
{
    int x_wall = grid->width / 8;
    int y_wall_start = 3*grid->height / 8;
    int y_wall_end =  5 * grid->height / 8;

    if(x == x_wall && y >= y_wall_start && y <= y_wall_end)
    {
        return 1;
    }
    return 0;
    /*
    int r = grid->height / 8;
    int cy = grid->height / 2;
    int cx = grid->width / 2;
    int dx = x - cx;
    int dy = y - cy;
    if(dx*dx + dy*dy < r*r)
    {
        return 1;
    }*/
        
    return 0;
}

void run_timestep(cell_grid* grid, double dt)
{
    double gamma = 1.4;
    cell_state_2d** new_cells = malloc(grid->height * sizeof(cell_state_2d*));
    for(int y = 0; y < grid->height; y++)
    {
        new_cells[y] = malloc(grid->width * sizeof(cell_state_2d));
        for(int x = 0; x < grid->width; x++)
        {
            cell_state_2d current = grid->cells[y][x];
            
            cell_state_2d top = {0};
            if(y == 0 || is_wall(x, y-1, grid))
            {   
                top = current;
                top.velocity.y = -current.velocity.y;
                top.momenta.y = -current.momenta.y;
            }
            else
            {
                top = grid->cells[y-1][x];
            }
            
            cell_state_2d right = {0};

            if(x == grid->width - 1)
            {
                right = current;
            }
            else if(is_wall(x+1, y, grid))
            {
                right = current;
                right.velocity.x = -current.velocity.x;
                right.momenta.x = -current.momenta.x;
            }
            else
            {
                right = grid->cells[y][x+1];
            }
            
            cell_state_2d bottom = {0};
            if(y == grid->height - 1 || is_wall(x, y+1, grid))
            {   
                bottom = current;
                bottom.velocity.y = -current.velocity.y;
                bottom.momenta.y = -current.momenta.y;
            }
            else
            {
                bottom = grid->cells[y+1][x];
            }
            cell_state_2d left = {0};
            
            double U_inlet = 1.0;
            double p_inlet = 2.0;
            double rho_inlet = 1.0;

            //enforce inlet velocity
            if(x == 0)
            {
                
                left.density = rho_inlet;
                left.velocity.x = U_inlet; 
                left.momenta.x = left.density * U_inlet;
                left.pressure = p_inlet;
                left.energy = (p_inlet / (gamma - 1)) + 0.5 * left.density * U_inlet * U_inlet;   
                
                //uncomment for outflow
                //left = current;
            }
            else if(is_wall(x-1, y, grid))
            {
                left = current;
                left.velocity.x = -current.velocity.x;
                left.momenta.x = -current.momenta.x;
            }
            else
            {
                left = grid->cells[y][x-1];
            }

            
            

            cell_state_2d new_state = compute(current, top, right, bottom, left, grid->dx, dt);
            new_cells[y][x] = new_state;
        }
    }
    for(int y = 0; y < grid->height; y++)    {
        for(int x = 0; x < grid->width; x++)
        {
            grid->cells[y][x] = new_cells[y][x];
        }
        free(new_cells[y]);
    }
    free(new_cells);
}


cell_grid * create_grid(int width, int height, double dx)
{
    double gamma = 1.4;
    cell_grid * grid = malloc(sizeof(cell_grid));
    grid->width = width;
    grid->height = height;
    grid->dx = dx;
    grid->cells = malloc(height * sizeof(cell_state_2d*));
    for(int y = 0; y < height; y++)
    {
        grid->cells[y] = malloc(width * sizeof(cell_state_2d));
        for(int x = 0; x < width; x++)
        {
            cell_state_2d state = {0};
            state.pressure = 2;
            state.density = 1.0;
            /*if(x < width/2)
            {
                state.density = 1.0;
                state.pressure = 1;
            }
            else
            {
                state.density = 0.125;
                state.pressure = 0.1;
            }*/
            
            state.velocity.x = 0;//((double)x / (double)width);
            state.velocity.y = 0.0;
            state.momenta.x = state.density * state.velocity.x;
            state.momenta.y = state.density * state.velocity.y;
            state.energy = (state.pressure / (gamma - 1)) + 0.5 * state.density * (state.velocity.x*state.velocity.x + state.velocity.y*state.velocity.y);
            grid->cells[y][x] = state;
        }
    }

    return grid;
}

void free_grid(cell_grid* grid)
{
    for(int y = 0; y < grid->height; y++)
    {
        free(grid->cells[y]);
    }
    free(grid->cells);
    free(grid);
}