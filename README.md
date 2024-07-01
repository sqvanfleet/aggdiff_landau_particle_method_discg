# Fully discrete energy-dissipative and conservative discrete gradient particle methods for a class of continuity equations
This repository includes the MATLAB code for each numerical example given in section 4 of 
*insert paper link*.  The goal is that any reader may reproduce the results shown for each example.  
## Contents 
There are five directories in this repository, four of them are aggregation-diffusion equation examples:
1. `particle_method_1D_heat_equation_discrete_gradient`
2. `particle_method_1D_porous_medium_discrete_gradient`
3. `particle_method_1D_linear_fokker_planck_discrete_gradient`
4. `particle_method_1D_non_local_fokker_planck_discrete_gradient`

The fifth directory Discrete_Gradient_Symmetric_reg, contains two subdirectories

1. `Discrete_Gradient_Symmetric_reg`
2. `Anisotropic_solution_with_Coulomb_potential`

Each of the six directories listed above contain two subdirecties:
1. `Plots`
2. `Data`

and the .m files:
1. `particle_method.m or particle_method_2d_parallel.m`
2. `gpsi_1d.m or gpsi_2d.m`
3. `right_hand_side.m or right_hand_side_parallel.m`
4. `plots.m`
5. `exact.m or exact_2d.m or non_bkw_initial_conditions.m`
6. `lgwt.m`

## MATLAB requirements

All of the files in this repository were tested on MATLAB version R2022b.  The examples coressponding to the Landau equation require the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html), but can be modified to be used in serial.  Each directory listed above contains a lgwt.m file that can be found 
at [MATLAB file exchange](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes).

## Project discription

1. `The particle.m or particle_method_2d_parallel.m` scripts begin by defining several parameters including:
    - `n_list` which is a matrix whose elements are how many particles you want to use.  The entries of         n_list correspond to different values of $M$ from the paper, making the total number of particles 
    $M^d$,    where $d$ is the number of spatial dimensions.  For example,
      ```matlab
      n_list = [60,70,80,90,100];
      ```
      means that the particle solution will be computed for $60^d, 70^d, 80^d, 90^d, 100^d$ particles.
    - 'Xmax' or 'Vmax' gives the computational domain which is the hypercube centered at the origin               
    $[\mbox{Xmax},\mbox{Xmax}]^d$ or $[\mbox{Vmax},\mbox{Vmax}]^d.
    - `dx` or `dv` is cell length so that the cell volume is $dx^d$ or $dv^d$.
    - The particle locations are initialized at the cell centers and for 1D examples we have
    ```matlab
    x = (-Xmax+dx/2):dx:(Xmax-dx/2);
    ```
    and for 2D examples we use
   ```matlab
   v = (-Vmax+dv/2):dv:(Vmax-dv/2);
   [vx,vy] = ndgrid(v);
   Vx = vx(:);
   Vy = vy(:);
   ```
   to create two arrays `Vx` and `Vy` which contain the corrdinates of the particles the $x$ and $y$ 
   direction.
   - `t0` is the initial time.
   - `f0` are the initial conditions which are computed using the `exact.m`, 
   `exact_2d.m`, or `non_bkw_initial_conditions.m` functions.
   - `w` are the particle weights and they are initialized using the midpoint    rule.
   - 
    
    



