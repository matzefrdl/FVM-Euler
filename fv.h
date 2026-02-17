#ifndef FV_H_
#define FV_H_

typedef struct vec2
{
    double x, y;
} vec2;

typedef struct cell_state
{
    double pressure;
    double velocity; // only one because of dim splitting
    double density;
} cell_state;

typedef struct cell_state_2d
{
    double pressure;
    vec2 velocity;
    vec2 momenta; // velocity * density
    double density;
    double energy;
} cell_state_2d;

#endif