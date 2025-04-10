# 1-D Advection-Diffusion Equation 

The 1-D advection-diffusion equation is a parabolic partial differential equation that is given by:

<p align="center">
  <img src="https://github.com/user-attachments/assets/c21ece9a-71d2-4ec7-8bae-0263048948dd" alt="Your Image">
</p>

where:

- **u(x, t)** is the scalar quantity being transported,
- **v** is the advection velocity,
- **D** is the diffusivity coefficient.

This equation describes the transport of a substance in a medium where both **advection** (transport due to bulk motion) and **diffusion** (spreading due to random motion) phenomoneon are happening simulateneously.

## Crank-Nicolson Method 

The **Crank-Nicolson method** is an implicit finite difference scheme that is unconditionally stable, meaning the numerical solution won’t blow up. However, stability alone doesn’t guarantee a physically meaningful solution, especially when advection is involved. The value of accuracy factors (r_diff and r_adv) can serve as guidelines for the user in defining the simulation parameters in [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt).   

The inital distribution of matter is modeled as a Gaussian pulse:

<p align="center">
  <img src="https://github.com/user-attachments/assets/02c713d8-1864-4da4-91ad-01988544d2eb" alt="Your Image">
</p>

Homogeneous Dirichlet boundary conditions are applied which enforce: 

<p align="center">
  <img src="https://github.com/user-attachments/assets/df167b36-6706-442a-b3be-79715a50d331" alt="Your Image">
</p>

Physically, this means that the concentration of matter that is being transported is **fixed at zero** at both ends of the domain for all times.

The spatial domain [0, L] and time domain [0, T] are divided into nx grid points and nx time steps respectively and then corresponding step sizes dx and dt are calculated. The discretization of the 1D advection-diffusion equation using the Crank-Nicolson method is then given as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/8e8b809c-101c-4024-bb8f-12fdaff727f6" alt="Your Image">
</p>

with the coefficients:

<p align="center">
  <img src="https://github.com/user-attachments/assets/3fa41534-fa45-4a0f-96aa-c81b47d84af9" alt="Your Image">
</p>

and the accuracy factors:

<p align="center">
  <img src="https://github.com/user-attachments/assets/275ed9bf-c7f4-4387-be9c-093fa9c84e65" alt="Your Image">
</p>

The equation can also be written in matrix form as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/7b26da2a-1573-4c20-86dd-e71ffc0724c8" alt="Your Image">
</p>

where A and B are tridiagnol matrices of size **nx X nx** with Dirichlet boundary conditions applied: 

 <p align="center">
  <img src="https://github.com/user-attachments/assets/05410a2d-d347-403f-9efa-79f0de9cc97c" alt="Your Image">
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/1bc51c0d-1915-48f6-abb6-44c7cf95c0e3" alt="Your Image">
</p>

## Analytical Solution 

The 1D Advection-Diffusion equation can be solved analytically by summing over a series of mirrored Gaussian pulses. The solution is expressed as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/44642763-139c-4010-af06-a4c1ae0df149" alt="Your Image">
</p>

Where:

**𝜎** is the standard deviation of the initial Gaussian pulse.
**𝐷** is the diffusivity coefficient.
**𝑣** is the advection velocity.
**𝐿** is the length of the spatial domain.
**𝑡** is the time.
**𝑥₀** is the initial center of the Gaussian pulse.

The term **2n𝐿** accounts for the mirrored reflections across the boundaries.

# Project Execution and Structure 

These are the steps to be followed to execute the program and get plotted results:

1. The user defines the simulation parameters in the [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt) file based on its given syntax or can use the pre-defined values. Note that the parameters: length, total_time, nx, nt, diffusivity and velocity are **required** by the program to run. The parameters x0, sigma, and num_reflections are optional — if not provided in the configuration file, the program assigns them sensible defaults: x0 is set to one-third of the domain length (length / 3), sigma to one-tenth of the domain length (length / 10), and num_reflections defaults to 5. These defaults are chosen to generate a well-behaved Gaussian pulse and ensure meaningful analytical solutions without requiring user input. The local paths for saving the numerical and analytical solutions must also be specified in the configuration file.
2. To start the solver, the user has to run [run_solver.py](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/run_solver.py) which imports the simulation parameters from [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt) by default using ConfigParser library. Note that it is possible to use a different file name than **configuration.txt** for the configuration file but the user has to specify the file name when launching the solver from the command line with the syntax ***"python run_solver.py custom_config.txt"***. If the user does not mention the name of the new file, the solver will read the simulation parameters from [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt) only. 
3. The analytical and numerical solutions are saved at the paths defined in the configuration file. The 1-D plot appears first comparing the numerical and analytical solutions over the same grid. The user can manually save the plot if desired. Closing that plot causes the 3D surface plot to appear showing the evolution of the numerical solution through space and time. Snapshots of the 3D solution can be saved from any desired perspective by adjusting the azimuthal and elevation angles to achieve the preferred viewpoint.

This is how I divided the project into blocks:

- In [functions.py](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/functions.py), I have built the functions needed for calculating discretization, calculating and checking accuracy factors as per guidelines, setting up inital gaussian pulse, creating matrices for C-N method, applying homogeneous dirichlet boundary conditions and solving the 1-D Advection-Diffusion equation numerically using the Crank-Nicolson method and analytically by summing over a series of mirrored Gaussian pulses.
- In [testing.py](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/testing.py), I have tested all the functions in functions.py to ensure that they work as expected.
- In [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt), I have pre-defined the simulation parameters needed for the program to run as well as the local paths where the numerical and analytical solution will be saved.
- In [run_solver.py](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/run_solver.py), this is the heart of the project that needs to be executed by the user. Here, the simulation parameters are extracted from the chosen configuration file, accuracy conditions are checked, and both numerical and analytical solutions are computed and saved. In the end, a 1D comparison plot and a 3D surface plot of the numerical solution are generated.
- In [plot.py](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/plot.py), there are the two functions that respectively plot the time evolution of the initial gaussian distribution as per the numerical and analytical solution for comparison and make the 3D surface plot of the numerical solution over space and time.

To show the results based on the current [configuration.txt](https://github.com/madhurpratap-s/1D-Advection-Diffusion-Solver/blob/main/configuration.txt):

1. **Plot comparing the time evolution of the initial distribution as per the numerical and analytical solution.**

Note that the advection phenomenon is apparent by the displacement of the inital gaussian distribution while the diffusion phenomenon can be observed in the form of broadening of the peaks over time along with reduction of the peak value.

![1D Advection-Diffusion Solution](./plots/1d_solution.png)

2. **3D surface plot of the numerical solution over space and time.**

![1D Advection-Diffusion Solution](./plots/3d_solution.png)
