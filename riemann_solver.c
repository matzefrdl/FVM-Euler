#include "riemann_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double f_X(double p_star, double rho_X, double p_X)
{
    double fX = 0;
    double gamma = 1.4;
    if(p_star > p_X)
    {
        double A = rho_X * (p_star * (gamma+1) + p_X * (gamma-1));
        fX = (p_star - p_X) * sqrt(2.0 / A);
    }
    else
    {
        double a_X = sqrt(gamma * p_X / rho_X);
        double pre = 2 * a_X / (gamma - 1.0);
        fX = pre * (pow(p_star/p_X, (gamma-1.0)/(2*gamma))-1);
    }

    return fX;
}

double pressure_root(double p_star, double p_r, double p_l, double u_r, double u_l, double rho_r, double rho_l)
{
    //fr + fl + ur - ul
    double fl = f_X(p_star, rho_l, p_l);
    double fr = f_X(p_star, rho_r, p_r);

    return fr + fl + u_r - u_l;
}

double find_p_star(double p_r, double p_l, double u_r, double u_l, double rho_r, double rho_l)
{
    double p_star = (p_l + p_r) / 2.0;

    double h = 1e-8;
    double iterations = 0;
     
    while(1)
    {
        double f = pressure_root(p_star, p_r, p_l, u_r, u_l, rho_r, rho_l);
        
        if(fabs(f) <= 1e-6)
        {
            break;
        }

        double fh = pressure_root(p_star+h, p_r, p_l, u_r, u_l, rho_r, rho_l);
        
        double df = (fh-f) / h;
        
        p_star = p_star - f / df;
        iterations++;
        if(iterations > 100)
        {
            //printf("Newton did not converge\n");
            
            return - 1234;
        }
    }

    return p_star;
}

riemann_sol solve(cell_state left, cell_state right)
{
    riemann_sol sol = {0};
    double gamma = 1.4;
    
    double p_star = find_p_star(right.pressure, left.pressure, right.velocity, left.velocity, right.density, left.density);
    
    if(p_star == -1234)
    {
        printf("-------------------\n");
        printf("U_l: %f, p_l: %f, rho_l: %f\n", left.velocity, left.pressure, left.density);
        printf("U_r: %f, p_r: %f, rho_r: %f\n", right.velocity, right.pressure, right.density);
        printf("-------------------\n");
        exit(1);
    }
    p_star = fmax(p_star, 0);
    sol.pressure_star = p_star;
    sol.velocity_star = right.velocity + f_X(p_star, right.density, right.pressure);

    sol.left_state = p_star > left.pressure ? 0 : 1; 
    sol.right_state = p_star > right.pressure ? 0 : 1; 
    if(sol.left_state) //rarefaction left
    {
        sol.rho_s_left = left.density * pow(p_star / left.pressure, 1.0 / gamma);
        
        double a_l = sqrt(gamma * left.pressure / left.density);
        double a_s_l = a_l * pow(p_star / left.pressure, (gamma - 1) / (2 * gamma));

        sol.SHL = left.velocity - a_l;
        sol.STL = sol.velocity_star - a_s_l;
    }
    else //shock left
    {
        double alpha = (gamma-1.0) / (gamma+1.0);
        sol.rho_s_left = left.density * ((p_star + alpha * left.pressure) / (left.pressure + alpha * p_star));
        
        double M_l = sqrt(- (p_star - left.pressure)/(1/sol.rho_s_left - 1/left.density));    
        sol.STL = left.velocity - M_l / left.density; // TODO: CHECK
        sol.SHL = left.velocity - M_l / left.density; // TODO: CHECK
    }
    if(sol.right_state) //rarefaction right
    {
        sol.rho_s_right = right.density * pow(p_star / right.pressure, 1.0 / gamma);
    
        double a_r = sqrt(gamma * right.pressure / right.density);
        double a_s_r = sqrt(gamma * p_star / sol.rho_s_right);

        sol.SHR = right.velocity + a_r;
        sol.STR = sol.velocity_star + a_s_r;
    
    }
    else //shock right
    {
        double alpha = (gamma-1.0) / (gamma+1.0);
        sol.rho_s_right = right.density * (p_star + alpha * right.pressure) / (right.pressure + alpha * p_star);
      
        double M_r = sqrt(- (right.pressure - p_star)/(1/right.density - 1/sol.rho_s_right));
        sol.SHR = right.velocity + M_r / right.density;
        sol.STR = right.velocity + M_r / right.density;
    }   

    return sol;
}


cell_state solve_new_state(riemann_sol sol, cell_state left, cell_state right)
{
    double gamma = 1.4;
    double SHX = 0;
    if(sol.velocity_star > 0)
    {  
        SHX = sol.SHL;
    
        if(sol.left_state) // rarefaction
        {
            if(SHX > 0) return left;
            if(sol.STL > 0) 
            {
                double a_L = sqrt(gamma * left.pressure / left.density);
                
                double vel = left.velocity * (gamma-1) / (gamma+1) + 2 / (gamma+1) * a_L;

                double p = left.pressure * pow(2 / (gamma + 1) + (gamma-1)/((gamma+1) * a_L)*left.velocity, 2.0 * gamma / (gamma - 1.0));
                
                double rho = left.density * pow(p / left.pressure, 1.0/gamma);

                cell_state state;
                state.pressure = p;
                state.velocity = vel;
                state.density = rho;

                return state;
            }
            else
            {
                cell_state state;
                state.pressure = sol.pressure_star;
                state.velocity = sol.velocity_star;
                state.density = sol.rho_s_left;
                return state;
            }
        }
        else //shock
        {
            if(0 < SHX) return left;
            if(SHX < 0)
            {
                cell_state state;
                state.pressure = sol.pressure_star;
                state.velocity = sol.velocity_star;
                state.density = sol.rho_s_left;
                return state;
            }
        }
    }
    else
    {  
        SHX = sol.SHR;
    
        if(sol.right_state)
        {
            if(SHX < 0)
            {
                return right;
            }
            if(sol.STR < 0)
            {
                double a_R = sqrt(gamma * right.pressure / right.density);
                
                double vel = right.velocity * (gamma-1) / (gamma+1) - 2 / (gamma+1) * a_R;

                double p = pow(vel * vel / gamma * right.density / pow(right.pressure, 1/gamma), gamma/(gamma-1));

                double rho = right.density * pow(p / right.pressure, 1.0/gamma);

                cell_state state;
                state.pressure = p;
                state.velocity = vel;
                state.density = rho;

                return state;
            }
            cell_state state;
            state.pressure = sol.pressure_star;
            state.velocity = sol.velocity_star;
            state.density = sol.rho_s_right;

            return state;
            
        }
        else //shock
        {
            if(0 > SHX) return right;
            if(SHX > 0)
            {
                cell_state state;
                state.pressure = sol.pressure_star;
                state.velocity = sol.velocity_star;
                state.density = sol.rho_s_right;
                return state;
            }
        }
    }
    cell_state state = {0};
    state.pressure = -1;
    return state;
}