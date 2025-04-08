# 1-D Advection-Diffusion Equation 

The 1-D advection-diffusion equation is a parabolic partial differential equation that is given by:

![image](https://github.com/user-attachments/assets/8737bf6d-6f5b-4a4b-8d07-a929c46526a1)


where:

- **u(x, t)** is the scalar quantity being transported,
- **v** is the advection velocity,
- **D** is the diffusivity coefficient.

This equation describes the transport of a substance in a medium where both **advection** (transport due to bulk motion) and **diffusion** (spreading due to random motion) phenomoneon are happening simulateneously.
___
## Crank-Nicolson Method 

The **Crank-Nicolson method** is an implicit finite difference scheme that is unconditionally stable, meaning the numerical solution won’t blow up. However, stability alone doesn’t guarantee a physically meaningful solution, especially when advection is involved. The value of accuracy factors (r_diff and r_adv) can serve as guidelines for the user in defining the simulation parameters in configuration.txt.   

The inital distribution of matter is modeled as a Gaussian pulse:

![image](https://github.com/user-attachments/assets/a8e2fc98-99b8-4755-b3f9-29ebb1393be3)

Homogeneous Dirichlet boundary conditions are applied which enforce: 

![image](https://github.com/user-attachments/assets/ce756ad6-e50c-4c46-b293-ff9f1b0ba30c)

Physically, this means that the concentration of matter that is being transported is **fixed at zero** at both ends of the domain for all times.

The spatial domain [0, L] and time domain [0, T] are divided into nx grid points and nx time steps respectively and then corresponding step sizes dx and dt are calculated. The discretization of the 1D advection-diffusion equation using the Crank-Nicolson method is then given as:

![image](https://github.com/user-attachments/assets/379c4acb-b35e-4cdb-aea9-92e766f95ad6)

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

![image](https://github.com/user-attachments/assets/c5137e32-98ea-4887-8376-8592691a880c)

where A and B are tridiagnol matrices of size **nx X nx** with Dirichlet boundary conditions applied: 

(Insert matrices A and B here later on)
___
## Analytical Solution 

The 1D Advection-Diffusion equation can be solved analytically by summing over a series of mirrored Gaussian pulses. The solution is expressed as:

![image](https://github.com/user-attachments/assets/1dfec78b-7f25-4907-8fb5-e28b2b7ecb72)

Where:

**𝜎** is the standard deviation of the initial Gaussian pulse.
**𝐷** is the diffusivity coefficient.
**𝑣** is the advection velocity.
**𝐿** is the length of the spatial domain.
**𝑡** is the time.
**𝑥₀** is the initial center of the Gaussian pulse.

The term **2m𝐿** accounts for the mirrored reflections across the boundaries.

# Project Structure 
