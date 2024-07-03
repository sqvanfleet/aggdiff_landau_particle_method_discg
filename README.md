# Fully discrete energy-dissipative and conservative discrete gradient particle methods for a class of continuity equations
This repository includes the MATLAB code for each numerical example given in section 4 of 
[this paper](https://arxiv.org/abs/2407.00533).  The goal is that any reader may reproduce the results shown for each example.  

## Collaborators 


## Contents 

There are five directories in this repository, four of them are aggregation-diffusion equation examples, which correspond to examples 4.1 - 4.4:
1. `particle_method_1D_heat_equation_discrete_gradient` (example 4.1)
2. `particle_method_1D_porous_medium_discrete_gradient` (example 4.2)
3. `particle_method_1D_linear_fokker_planck_discrete_gradient` (example 4.3)
4. `particle_method_1D_non_local_fokker_planck_discrete_gradient` (example 4.4)

The fifth directory Discrete_Gradient_Symmetric_reg, contains two subdirectories, which correspond to examples 4.5 and 4.6:

1. `Discrete_Gradient_Symmetric_reg` (example 4.5)
2. `Anisotropic_solution_with_Coulomb_potential` (example 4.6)

Each of the six directories listed above contain two subdirecties:
1. `Plots` this directory stores all the plots created by the `plots.m` file 
2. `Data` this directory stores all of the data saved from the `particle_method.m` or `particle_method_2d_parallel.m` file

and the .m files:
1. `particle_method.m or particle_method_2d_parallel.m` the main script that preforms the computations and saves data to the `Data` directory
    that is then analyzed with the `plots.m` file
3. `gpsi_1d.m or gpsi_2d.m` functions for the mollifier funtion $\varphi_{\varepsilon}(\boldsymbol{x})$
4. `right_hand_side.m or right_hand_side_parallel.m` a matlab function that computes the right hand side of the fixed point
   iteration method that results from the discrete gradient discretization.  
5. `plots.m` takes data from `Data` and produces pdf plots that are saved in `Plots`
6. `exact.m or exact_2d.m or non_bkw_initial_conditions.m` functions for the exact solution for the given example or for example 4.6
   a function for the initial conditions
8. `lgwt.m` a function that returns the wights and nodes required for the Gauss-Legendre quadrature.  This can be downloaded at
   at [MATLAB file exchange](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes).

## MATLAB requirements

All of the files in this repository were tested on MATLAB version R2022b.  The examples coressponding to the Landau equation require the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html), but can be modified to be used in serial. 

## Project discription

### `particle_method.m` or `particle_method_2d_parallel.m` scripts

#### Parameters, initialization, and mesh generation

- The `particle_method.m` or `particle_method_2d_parallel.m` scripts begin by defining several parameters including:
- `n_list` which is a matrix whose elements are how many particles you want to calculate the particle method solution with.  The entries of this array are chosen by the user.  The entries of `n_list` correspond to different values of $M$ from the paper, making the total number of particles $M^d$, where $d$ is the number of spatial dimensions.  The outer most loop in `particle_method.m` or `particle_method_2d_parallel.m` uses the loop variable `alpha` and loops over all values in `n_list`.  For each iteration of this loop the number of particles is
```matlab
n = n_list(alpha);
```
- For the 2D examples,
```matlab
Np = n^2
```
is the number of particles.  

- Landau equation examples have two additional parameters `gamma` and `C_gamma`, which are parameters representing the collision kernal 
$$a_{ij}(\boldsymbol{x})=C|\boldsymbol{x}|^{\gamma}(|\boldsymbol{x}|^2\delta_{ij}-x_ix_j)$$

- `Xmax` or `Vmax` gives the computational domain which is the hypercube centered at the origin               
$[\mbox{Xmax},\mbox{Xmax}]^d$ or $[\mbox{Vmax},\mbox{Vmax}]^d$. These variables are chosen by the user.
- `dx` or `dv` is cell length so that the cell volume is $dx^d$ or $dv^d$, and is computed with
```matlab
dx = 2*Xmax/n
```
```matlab
dv = 2*Vmax/n
```
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
- `t0` is the initial time and chosen by the user.
- `f0` or `f_initial` are the initial conditions which are computed using the `exact.m`, `exact_2d.m`, or `non_bkw_initial_conditions.m` functions.  
- `w` are the particle weights and they are initialized using the midpoint rule.  In the 2D examples `W = w(:)` flattens the `w` array so that `W` is the same size as the `Vx` and `Vy` arrays. 
- To create the "blob solution" (equation 4.1), the reconstruction mesh is created.  For 1D exampels
```matlab
Nr = n;
dxr = 2*Xmax/n; 
xr = (-Xmax+dxr/2):dxr:(Xmax-dxr/2); 
```
   
and for 2D examples
```matlab
Nr = n;
dvr = 2*Vmax/Nr;
vr = (-Vmax+dvr/2):dvr:(Vmax-dvr/2); 
[vrx,vry] = ndgrid(vr);
Vrx = vrx(:);
Vry = vry(:);
```
Creating this reference mesh is required to compute the errors or plot the blob solution.
- Once the meshsize `dx` or `dv` is choses the regularization parameter $\epsilon$ is chosen to be
```matlab
epsilon = 4*(0.4*(dv)^0.99)^2;
```
- The initial "reconstructed or blob solution" can then be calculated using the formula $$\sum^N_{p=1} w_p\varphi_{\varepsilon}(\boldsymbol{x} - \boldsymbol{x}_p).$$
###   The `psi_1d.m` and `psi_2d.m` functions
- The arguments for the `psi_1d` are `x` and `eps` and for `psi_2d` they are `vx`, `vy`, and `eps`.
- These funtions return an array the same size as `x` or `vx` and `vy`.  Each entry of the retuned array is the mollifier function
  $$\varphi_{\varepsilon}(\boldsymbol{x}) = \frac{1}{2 \pi \varepsilon}\exp{\left(\frac{-|\boldsymbol{x}^2|}{2 \varepsilon}\right)},$$
  evauluated at the corresponding element of `x` or `vx` and `vy` with the parameter `eps`.
  
```matlab
f = zeros(1,Nr);
for i = 1:Nr
    f(i) = sum(w.*psi_1d(xr(i)-x,epsilon));
end
```
and in 2D
```matlab
parfor i = 1:Nr
    for j = 1:Nr
        f(i,j) = sum(W.*psi_2d(vrx(i,j)-Vx,vry(i,j)-Vy,epsilon));
    end
end
```
- `tmax` is the final time and is chosen by the user.
- `dt` is the time step and is chosen by the user.
- `Nt` is the number of time steps
  ```matlab
  Nt = round((tmax-t0)/dt);
  ```
- `error_list` is an array that stores several quantities at each time step such as errors, mass, momentum, kinetic energy,
   energy, the dissipation term, fisher energy, and how many fixed point iterations were required for the relative error to
   be less than the given tolerance.
- `max_iter` is the maximum amount of iterations the fixed point iteration can preform and is chosen by the user.
- `tol` is the tolerence used in the fixed point iteration and is chosen by the user.

#### For loop in time

- The for loop in time uses `nt` as the loop variable and goes from 1 to `Nt`
  
- `time` represents the current time
```matlab
time = t0+dt*nt
```
- Next, a forward Euler step is calculated and used for the initial guess in the fixed point iteration

    



