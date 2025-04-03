# 1-D Advection-Diffusion Equation Solver

This repository contains a Python implementation of a **1-D Advection-Diffusion equation solver** using the **Crank-Nicolson method**. The solver computes the numerical solution and compares it with the analytical solution obtained using the **Gaussian Mirror Method**.

## 📖 1-D Advection-Diffusion Equation

The 1-D advection-diffusion equation is given by:

\[
\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = D \frac{\partial^2 u}{\partial x^2}
\]

where:

- \( u(x,t) \) is the scalar quantity being transported,
- \( v \) is the advection velocity,
- \( D \) is the diffusion coefficient.

This equation describes the transport of a substance in a medium with both **advection** (transport due to bulk motion) and **diffusion** (spreading due to random motion).

## ⚙️ Crank-Nicolson Method

The **Crank-Nicolson method** is an implicit finite difference scheme that is unconditionally stable, meaning the numerical solution won’t blow up. However, stability alone doesn’t guarantee a physically meaningful solution — especially when advection is involved. This is addressed in the code through calculation of accuracy factors whose values can guide the user in defining the configuration file.  

The Crank-Nicolson method uses a time-centered difference to discretize the equation:

\[
\frac{u_i^{n+1} - u_i^n}{\Delta t} + v \frac{u_{i+1}^{n+1} - u_{i-1}^{n+1} + u_{i+1}^{n} - u_{i-1}^{n}}{4 \Delta x} = D \frac{u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1} + u_{i+1}^{n} - 2u_i^{n} + u_{i-1}^{n}}{2 \Delta x^2}
\]

This leads to a **tridiagonal system of equations** that can be solved efficiently using the Thomas algorithm.

## 📊 Analytical Solution Using the Gaussian Mirror Method

The **Gaussian Mirror Method** provides an analytical solution to the **1D Advection-Diffusion Equation** for a finite domain by considering mirrored image sources.

### 🔢 Mathematical Formulation

For an initial Gaussian pulse centered at \( x_0 \) with standard deviation \( \sigma \), the free-space solution for the advection-diffusion equation is:

\[
u(x,t) = \frac{\sigma}{\sqrt{\sigma^2 + 2Dt}} \exp\left( - \frac{(x - (x_0 + vt))^2}{2 (\sigma^2 + 2Dt)} \right)
\]

However, in a **finite domain** with **Dirichlet boundary conditions**, the **Gaussian Mirror Method** accounts for boundaries by summing over reflected images of the Gaussian pulse:

\[
u(x,t) = \sum_{n=-N}^{N} \frac{\sigma}{\sqrt{\sigma^2 + 2Dt}} \exp\left( - \frac{(x - (x_0 + vt) - 2nL)^2}{2 (\sigma^2 + 2Dt)} \right)
\]

where:

- \( L \) is the domain length,
- \( v \) is the advection velocity,
- \( D \) is the diffusion coefficient,
- \( n \) represents the mirrored reflections,
- \( N \) is the number of reflections considered.

As \( N \to \infty \), this summation enforces the correct boundary conditions while preserving the Gaussian profile's evolution.
