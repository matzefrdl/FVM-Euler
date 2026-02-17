# Finite Volume Solver for Euler Equations

This code solves the 2D Euler equations. It uses [raylib](https://www.raylib.com/) for visualization.  
It currently supports only square cells as finite volumes.

## Simulation Setup
By default, a stream with a predefined velocity enters from the left.  
On the right, there is an outflow. On the top and bottom, walls are placed.  
Obstacles can be placed inside the domain.

<p align="center">
<img src="https://github.com/matzefrdl/FVM-Euler/blob/main/setup.png" alt="Setup for the simulation" width="800"/>  
</p>

## Solver
The code includes a Riemann solver for the 1D Euler equations. The 2D Euler equations are solved using dimensional splitting.  
The transverse velocity is handled as a passive scalar. For each cell, four Riemann problems are solved.

## Build/Run
**Prerequisites:** gcc, make, raylib

```console
[user@machine folder]$ make
gcc -g -Wall -Werror -o main discrete.c fv.c main.c riemann_solver.c -lm -lraylib 
[user@machine folder]$ ./main
```
## Result
A window will open with the dimensions of the mesh. By default, pressure will be shown. Pressing the following keys switches to a different field:
| Key  | Field shown |
| --- | --- |
| P  | Pressure |
| V  | Speed |
| D  | density |

The solver suffers from numerical dissipation and splitting errors.
As an example, the setup described above at T = 2:
<p align="center">
<img src="https://github.com/matzefrdl/FVM-Euler/blob/main/T2.png" alt="Setup described above with at T = 2" width="800"/>  
</p>
