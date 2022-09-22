This code runs a quasi-potential flow model for an impacting droplet on a bath of the same fluid.
The process is user-friendly, and the process is as follows:

1) input physical parameters and simulation settings into the driver_3B.m file. 
2) run the driver_3B.m file
3) This file loops over user-inputted impact velocities are repeatedly calls the bounce_alventosa_bessel.m function, which houses the solver, as well as the
  pre-computation of pressure distributions. If the pressure functions .mat files are not in the current working directory, the code will make them. 
  the bounce_alventosa_bessel.m function initializes the solver, and propogates through time until the user-input end time is reached. It then computes output impact parameters like
  coefficient of restitution, contact time, and maximum penetration depth. The intersections.m file is called within the solver to determine the instantaneous
  size of the pressure distribution, and the besselzero.m file determines zeros of bessel functions.
4) The driver code then calls the plot_drop_waves_3.m file, which creates 2 videos, one of the bath mode amplitudes as a function of time 
  (can toggle to plot drop mode amplitudes as a functiobn of time as well) and an  aesthetic video of the impact, either in 2D or 3D 
  (can toggle 2D or 3D mode in function input)
 5) that's about it!
