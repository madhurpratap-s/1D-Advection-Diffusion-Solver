# 1-D Advection-Diffusion Equation 

The 1-D advection-diffusion equation is a parabolic partial differential equation that is given by:

![image](https://github.com/user-attachments/assets/330bd080-156e-4205-a7df-fbd5ab97a4fa)

where:

- **u(x, t)** is the scalar quantity being transported,
- **v** is the advection velocity,
- **D** is the diffusivity coefficient.

This equation describes the transport of a substance in a medium where both **advection** (transport due to bulk motion) and **diffusion** (spreading due to random motion) phenomoneon are happening simulateneously.
___
## Crank-Nicolson Method 

The **Crank-Nicolson method** is an implicit finite difference scheme that is unconditionally stable, meaning the numerical solution won’t blow up. However, stability alone doesn’t guarantee a physically meaningful solution, especially when advection is involved. The value of accuracy factors (r_diff and r_adv) can serve as guidelines for the user in defining the simulation parameters in configuration.txt.   

The inital distribution of matter is modeled as a Gaussian pulse:

![image](https://github.com/user-attachments/assets/ade100dd-2d15-4e69-b892-950fa0e68f11)

Homogeneous Dirichlet boundary conditions are applied which enforce: 

![image](https://github.com/user-attachments/assets/daebeb1b-8ee1-4f91-8409-42547aa3431f)

Physically, this means that the concentration of matter that is being transported is **fixed at zero** at both ends of the domain for all times.

The spatial domain [0, L] and time domain [0, T] are divided into nx grid points and nx time steps respectively and then corresponding step sizes dx and dt are calculated. The discretization of the 1D advection-diffusion equation using the Crank-Nicolson method is then given as:

![image](https://github.com/user-attachments/assets/c9baab12-484b-4965-a54f-ae2d47eefdf2)

with the coefficients:

a_c = 1 + r_diff  
a_l = -0.5 r_diff - 0.25 r_adv  
a_u = -0.5 r_diff + 0.25 r_adv  
b_c = 1 - r_diff  
b_l = 0.5 r_diff + 0.25 r_adv  
b_u = 0.5 r_diff - 0.25 r_adv

and the accuracy factors:

r_diff = D * dt / dx²  
r_adv = v * dt / dx

The equation can also be written in matrix form as:

![image](https://github.com/user-attachments/assets/f1363ff1-c727-4d60-9bf5-bebf7d00e6f6)

where A and B are tridiagnol matrices of size **nx X nx** with Dirichlet boundary conditions applied: 

(Insert matrices A and B here later on)
___
## Analytical Solution 

The 1D Advection-Diffusion equation can be solved analytically by summing over a series of mirrored Gaussian pulses. The solution is expressed as:

![image](https://github.com/user-attachments/assets/f5128070-7d50-4f4b-8a79-3a8c7a548c5c)

Where:

**𝜎** is the standard deviation of the initial Gaussian pulse.
**𝐷** is the diffusivity coefficient.
**𝑣** is the advection velocity.
**𝐿** is the length of the spatial domain.
**𝑡** is the time.
**𝑥₀** is the initial center of the Gaussian pulse.

The term **2m𝐿** accounts for the mirrored reflections across the boundaries.

# Project Structure 
