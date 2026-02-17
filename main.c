#include "discrete.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include "raylib.h"

int main()
{
    cell_grid * grid = NULL;

    grid = create_grid(200, 200, 0.01);

    const int screenWidth = grid->width * 4;
    const int screenHeight = grid->height * 4;

    InitWindow(screenWidth, screenHeight, "Finite Volume Method");

    //SetTargetFPS(60);
    double T = 0;
    double dt = 0.001;
    int state_to_view = 0;

    int running = 0;
    while (!WindowShouldClose())
    {    
        T += dt;
        if(running)
        {
            printf("Time: %f\n", T);
            run_timestep(grid, dt);
        }
        if(IsKeyPressed(KEY_SPACE))
        {
            running = 1;
        }
        if(IsKeyPressed(KEY_P))
        {
            state_to_view = 0;
        }
        else if(IsKeyPressed(KEY_V))
        {
            state_to_view = 1;
        }
        else if(IsKeyPressed(KEY_D))
        {
            state_to_view = 2;
        }
        BeginDrawing();

            ClearBackground(BLACK);

            for(int y = 0; y < grid->height; y++)
            {
                for(int x = 0; x < grid->width; x++)
                {
                    cell_state_2d state = grid->cells[y][x];
                    
                    double var = 0;
                    if(state_to_view == 0)
                    {
                        var = state.pressure / 10.0;
                    }
                    else if(state_to_view == 1)
                    {
                        double speed = sqrt(state.velocity.x*state.velocity.x + state.velocity.y*state.velocity.y);
                        var = speed;
                    }
                    else if(state_to_view == 2)
                    {
                        var = state.density / 2.0;
                    }
                    int color_value = (int)(var * 255);
                    if(color_value > 255) color_value = 255;
                    Color color = {color_value, color_value, color_value, 255};
                    DrawRectangle(x*4, y*4, 4, 4, color);
                }
            }
            
        EndDrawing();
    }

    CloseWindow();

    free_grid(grid);

    return 0;
}
/* 1D SOD Shock tube problem
int main()
{

    double dx = 0.01;
    int N_cells = 100;
    double aspect_ratio = 1;
    const int screenWidth = N_cells * 4;
    const int screenHeight = N_cells * aspect_ratio * 4;

    InitWindow(screenWidth, screenHeight, "FPS");

    cell_state * grid = malloc(N_cells * sizeof(cell_state));
    double * momenta = malloc(N_cells * sizeof(double));
    double * energys = malloc(N_cells * sizeof(double));
    
    double gamma = 1.4;
    for(int i = 0; i < N_cells; i++)
    {
        cell_state state = {0};
        if(i < N_cells/2)
        {
            state.density = 1.0;
            state.pressure = 1.0;
        }
        else        {
            state.density = 0.125;
            state.pressure = 0.1;
        }
        state.velocity = 0.0;
        
        energys[i] = state.pressure / (gamma - 1) + 0.5 * state.density * state.velocity * state.velocity;
        momenta[i] = state.density * state.velocity;
        grid[i] = state;
    }

    SetTargetFPS(60);
    double T = 0;
    double dt = 0.001;
    while (!WindowShouldClose())
    {   
        if(T < 0.2) {
            T += dt;
        printf("Time: %f\n", T);
        

        cell_state * new_grid = malloc(N_cells * sizeof(cell_state));
        for(int i = 0; i < N_cells; i++)
        {
            //printf("%i\n", i);

            cell_state left = {0};
            if(i > 0)
            {
                left = grid[i-1];
            }
            else
            {
                left = grid[i];
            }
            cell_state right = i < N_cells - 1 ? grid[i+1] : grid[i];

            new_grid[i] = grid[i];
            riemann_sol sol_l = solve(left, grid[i]);
            riemann_sol sol_r = solve(grid[i], right);
            cell_state sol_r_state = solve_new_state(sol_r, grid[i],right);
            cell_state sol_l_state = solve_new_state(sol_l, left, grid[i]);

            double drhodt_l = -sol_l_state.density * sol_l_state.velocity;
            double dudt_l = -sol_l_state.density * sol_l_state.velocity * sol_l_state.velocity - sol_l_state.pressure;
            double u2_l = sol_l_state.velocity * sol_l_state.velocity;
            double E_l = 0.5 * sol_l_state.density * (u2_l) + sol_l_state.pressure / (gamma-1);
            double dEdt_l = -sol_l_state.velocity * (E_l + sol_l_state.pressure);
            
            double drhodt_r = -sol_r_state.density * sol_r_state.velocity;
            double dudt_r = -sol_r_state.density * sol_r_state.velocity * sol_r_state.velocity - sol_r_state.pressure;
            double u2_r = sol_r_state.velocity * sol_r_state.velocity;
            double E_r = 0.5 * sol_r_state.density * (u2_r) + sol_r_state.pressure / (gamma-1);
            double dEdt_r = -sol_r_state.velocity * (E_r + sol_r_state.pressure);

            new_grid[i].density += (drhodt_r - drhodt_l) * dt / dx;
            energys[i] += (dEdt_r - dEdt_l) * dt / dx;
            momenta[i] += (dudt_r - dudt_l) * dt / dx;
            new_grid[i].velocity = momenta[i] / new_grid[i].density;
            new_grid[i].pressure = (gamma - 1) * (energys[i] - 0.5 * new_grid[i].density * new_grid[i].velocity * new_grid[i].velocity);
            if(new_grid[i].pressure < 0)
            {
                printf("............... ERROR Report ..........\n");
                printf("Negative pressure at cell %i: %f\n", i, new_grid[i].pressure);
                printf("Density: %f, velocity: %f, energy: %f\n", new_grid[i].density, new_grid[i].velocity, energys[i]);
                printf("Left state: density: %f, velocity: %f, pressure: %f\n", left.density, left.velocity, left.pressure);
                printf("Right state: density: %f, velocity: %f, pressure: %f\n", right.density, right.velocity, right.pressure);
                exit(1);
            }
        }
        
        for(int i = 0; i < N_cells; i++)
        {
            grid[i] = new_grid[i];
        }
        free(new_grid);
    }

        BeginDrawing();

            ClearBackground(BLACK);
            
            //plot velocity
            for(int i = 0; i < N_cells; i++)
            {
                cell_state state = grid[i];
                double speed = state.density;

                double scaled_speed = 1-speed / 1.0;

                double y = scaled_speed * (screenHeight - 20) + 10;
                DrawCircle(i*4 + 2, (int)y, 2, WHITE);
            }

        EndDrawing();
    }
    free(grid);
    CloseWindow();
   


    return 0;
}
*/