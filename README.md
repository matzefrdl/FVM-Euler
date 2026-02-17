# Finite Volume Solver for Euler equations
This code solves the 2D euler equations. It uses [raylib](https://www.raylib.com/) for visualization.
It only supports squares as finite volumes. 

## Simulation Setup
By default on the left enters a stream with a predefined velocity.
On the right there is an outflow. On the top and bottom a wall is placed.
Obstacles can be placed into the volume.

## Solver
The code includes an riemann solver of the 1D euler equation. The 2D euler equation are solved with dimensional splitting. The transverse velocity is handled as a passive scalar. For each cell four riemann problems are solved.

## Build/Run
Prerequisites: gcc, make, raylib
```console
[user@machine folder]$ make
gcc -g -Wall -Werror -o main discrete.c fv.c main.c riemann_solver.c -lm -lraylib 
[user@machine folder]$ ./main
```
## Result
A window will open with the dimensions of the mesh. By default pressure will be shown. Pressing the following keys switches to a different field.
| Key  | Field shown |
| --- | --- |
| P  | Pressure |
| V  | Speed |
| D  | density |

As an example the setup describe above at T = ...:
<p align="center">
<img src="https://github.com/matzefrdl/BNN-with-HMC/blob/main/moons_hmc.png" alt="Make moons from sklearn classified with a BNN" width="500"/>  
</p>
