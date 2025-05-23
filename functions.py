# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:52:20 2025

@author: madhu
"""

import warnings
import numpy as np

def calculate_discretization(L, T, nx, nt):
    """
    Calculates spatial and temporal step sizes.
    
    This function computes the step sizes based on the provided domain length (L), total simulation time (T), 
    number of spatial steps (nx), and number of time steps (nt).

    Args:
        L (float): Length of the domain.
        T (float): Total simulation time.
        nx (int): Number of spatial steps.
        nt (int): Number of time steps.

    Returns:
        tuple: A tuple containing:
            - dx (float): Spatial step size.
            - dt (float): Temporal step size.
    
    Raises:
        ValueError: If nx or nt are less than 2, as they must represent valid discretization parameters.
    """
    if nx < 2 or nt < 2:
       raise ValueError("nx and nt must be greater than or equal to 2.")
    
    dx = L / (nx - 1)
    dt = T / (nt - 1)
    
    return dx, dt

def calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity):
    """
    Calculates and checks if accuracy factors are consistent with guidelines.

    Args:
        L (float): Length of the domain.
        T (float): Total simulation time.
        nx (int): Number of spatial steps.
        nt (int): Number of time steps.
        D (float): Diffusivity coefficient of the medium.
        velocity (float): Advection velocity.

    Returns:
        tuple: A tuple containing:
            - r_diff (float): Diffusion accuracy factor (D * dt / dx²).
            - r_adv (float): Advection accuracy factor (velocity * dt / dx).
    
    Warns:
        UserWarning: If the accuracy conditions are not met (if r_diff > 0.5 or r_adv > 1.0).
    """
    dx, dt = calculate_discretization(L, T, nx, nt)

    r_diff = D * dt / dx**2
    r_adv = velocity * dt / dx
    
    if r_diff > 0.5:
        warnings.warn(f"Potential accuracy issue: r_diff={r_diff}. Recommended r_diff < 0.5 for accuracy.", UserWarning)
    if r_adv > 1.0:
        warnings.warn(f"Potential accuracy issue: r_adv={r_adv}. Recommended r_adv < 1.0 for accuracy.", UserWarning)
    
    return r_diff, r_adv
    
def setup_gaussian_pulse(L, nx, x0=None, sigma=None):
    """
    Generates an initial Gaussian pulse distribution on a 1D grid.

    Args:
        L (float): Length of the domain.
        nx (int): Number of spatial steps.
        x0 (float, optional): Center of the Gaussian pulse. Defaults to L/3.
        sigma (float, optional): Standard deviation. Defaults to L/10.

    Returns:
        numpy.ndarray: Initial concentration profile.    
    """
    if x0 is None:
        x0 = L / 3  
    if sigma is None:
        sigma = L / 10  
        
    dx, dt = calculate_discretization(L, T=1.0, nx=nx, nt=5) # Use dummy T and any nt > 2
    x = np.linspace(0, (nx - 1) * dx, nx)
    
    return np.exp(-0.5 * ((x - x0) / sigma)**2)

def create_matrices(nx, r_diff, r_adv):
    """
    Calculates and returns the matrices A and B for the Crank-Nicolson method.
 
    Args:
        nx (int): Number of spatial grid points.
        r_diff (float): Diffusion number (D * dt / dx²).
        r_adv (float): Advection number (v * dt / dx).

    Returns:
        tuple: A tuple containing:
            - A (numpy.ndarray): Left-hand side matrix for the C-N time-stepping scheme.
            - B (numpy.ndarray): Right-hand side matrix for the C-N time-stepping scheme. 
    """
    # Main diagonals
    A_main = np.ones(nx) * (1 + r_diff)
    B_main = np.ones(nx) * (1 - r_diff)

    # Upper and lower diagonals 
    A_upper = np.ones(nx - 1) * (-0.5 * r_diff + 0.25 * r_adv)  
    A_lower = np.ones(nx - 1) * (-0.5 * r_diff - 0.25 * r_adv)  
    B_upper = np.ones(nx - 1) * (0.5 * r_diff - 0.25 * r_adv) 
    B_lower = np.ones(nx - 1) * (0.5 * r_diff + 0.25 * r_adv)  

    # Construct matrices
    A = np.diag(A_main) + np.diag(A_upper, k=1) + np.diag(A_lower, k=-1)
    B = np.diag(B_main) + np.diag(B_upper, k=1) + np.diag(B_lower, k=-1)

    return A, B

def apply_boundary_conditions(A, B): 
    """
    Applies homogeneous Dirichlet boundary conditions to matrices A and B.

    Args:
        A (numpy.ndarray): Left-hand side matrix for the Crank-Nicolson method.
        B (numpy.ndarray): Right-hand side matrix for the Crank-Nicolson method.

    Returns:
        tuple: A tuple containing:
            - A (numpy.ndarray): Modified left-hand side matrix with boundary conditions applied.
            - B (numpy.ndarray): Modified right-hand side matrix with boundary conditions applied.
    """
    A[0, :] = 0
    A[0, 0] = 1
    A[-1, :] = 0
    A[-1, -1] = 1

    B[0, :] = 0
    B[0, 0] = 1
    B[-1, :] = 0
    B[-1, -1] = 1

    return A, B

def solve_advection_diffusion_CN(L, T, nx, nt, D, velocity, x0, sigma):
    """
    Solves the 1D advection-diffusion equation using the Crank-Nicolson method.

    Args:
        L (float): Length of the domain.
        T (float): Total simulation time.
        nx (int): Number of spatial steps.
        nt (int): Number of time steps.
        D (float): Diffusivity coefficient of the medium.
        velocity (float): Advection velocity.
        x0 (float, optional): Center of the initial Gaussian pulse. Defaults to L/3.
        sigma (float, optional): Std. Dev. of the initial Gaussian Pulse. Defaults to L/10.

    Returns:
        tuple: A tuple containing:
            - x (numpy.ndarray): Spatial grid points of length `nx`.
            - u (numpy.ndarray): Concentration profile of shape [nx, nt] over time.   
    """
    r_diff, r_adv = calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)
    x = np.linspace(0, L, nx)
    u = np.zeros((nx, nt))
    u[:, 0] = setup_gaussian_pulse(L, nx, x0 = x0, sigma = sigma)

    # Apply Dirichlet boundary condition to solution array
    u[0, :] = 0
    u[-1, :] = 0

    A, B = create_matrices(nx, r_diff, r_adv)
    A, B = apply_boundary_conditions(A, B)

    # Time-stepping loop
    for n in range(1, nt):
        b = B @ u[:, n - 1]
        b[0] = 0
        b[-1] = 0
        u[:, n] = np.linalg.solve(A, b)

    return x, u

def solve_advection_diffusion_analytical(L, T, nx, nt, D, velocity, x0, sigma, num_reflections = 5):
    """
    Solves the 1D advection-diffusion equation using the gaussian mirror method based on the analytical solution.

    While the method leverages the analytical solution, it includes a numerical component 
    necessary for practical implementation. This component allows for plotting and comparing 
    the results with numerical methods like the Crank-Nicolson method.
    
    The solution is given by summing over mirrored Gaussians:

        u(x,t) = Σ (sigma / sqrt(sigma^2 + 2*D*t)) * exp( - (x - (x0 + velocity*t) - 2nL)^2 / (2*(sigma^2+2*D*t)) )

    Args:
        L (float): Length of the spatial domain.
        T (float): Total simulation time.
        nx (int): Number of spatial points.
        nt (int): Number of time steps.
        D (float): Diffusivity coefficient.
        velocity (float): Advection velocity.
        x0 (float, optional): Center of the initial Gaussian pulse. Defaults to L/3.
        sigma (float, optional): Std. Dev. of the initial Gaussian Pulse. Defaults to L/10.
        num_reflections (int, optional): Number of mirrored sources to consider. Defaults to 5.
    
    Returns:
        tuple: A tuple containing:
            - x (numpy.ndarray): Spatial grid points.
            - u (numpy.ndarray): 2D solution array (shape: (nx, nt)), with each column representing a time step.
    Raises:
        ValueError: If num_reflections is less than 1.
    """
    if num_reflections < 1:
        raise ValueError("num_reflections must be a positive integer (>= 1).")
        
    dx, dt = calculate_discretization(L, T, nx, nt)
    x = np.linspace(0, (nx - 1) * dx, nx)

    # Set default Gaussian pulse parameters
    if x0 is None:
        x0 = L / 3
    if sigma is None:
        sigma = L / 10

    u = np.zeros((nx, nt))

    # Time-stepping loop
    for n in range(nt):
        t = n * dt
        denom = sigma**2 + 2 * D * t

        for m in range(-num_reflections, num_reflections + 1):
            x_shift = x0 + velocity * t + 2 * m * L  # Mirrored image shifts
            u[:, n] += (sigma / np.sqrt(denom)) * np.exp(-((x - x_shift) ** 2) / (2 * denom))

    return x, u
