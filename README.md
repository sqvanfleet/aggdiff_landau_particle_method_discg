# Fully discrete energy-dissipative and conservative discrete gradient particle methods for a class of continuity equations
This repository includes the MATLAB code for each numerical example given in section 4 of 
*insert paper link*.  The goal is that any reader may reproduce the results shown for each example.  
## Contents 
There are five directories in this repository, four of them are aggregation-diffusion equation examples:
1. particle_method_1D_heat_equation_discrete_gradient
2. particle_method_1D_porous_medium_discrete_gradient
3. particle_method_1D_linear_fokker_planck_discrete_gradient
4. particle_method_1D_non_local_fokker_planck_discrete_gradient <br>

The fifth directory Discrete_Gradient_Symmetric_reg, contains two subdirectories

1. Discrete_Gradient_Symmetric_reg
2. Anisotropic_solution_with_Coulomb_potential

Each of the six directories listed above contain two subdirecties:
1. Plots
2. Data

and the .m files:
1. particle_method.m or particle_method_2d_parallel.m
2. gpsi_1d.m or gpsi_2d.m
3. right_hand_side.m or right_hand_side_parallel.m
4. plots.m
5. exact.m or exact_2d.m or non_bkw_initial_conditions.m
6. lgwt.m
