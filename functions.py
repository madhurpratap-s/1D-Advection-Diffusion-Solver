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
    
    This function computes the step sizes based on the provided domain length, total simulation time, 
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

def calculate_accuracy_factors(L, T, nx, nt, D, velocity):
    """
    Calculates the accuracy factors for the Crank-Nicolson method applied to the advection-diffusion equation.

    Although the Crank-Nicolson method is unconditionally stable for the advection-diffusion equation,
    poor choices of spatial and temporal discretization can lead to numerical inaccuracies such as
    oscillations or excessive numerical diffusion. These dimensionless numbers help assess the 
    suitability of discretization choices for maintaining accuracy in the simulation.

    Args:
        L (float): Length of the domain.
        T (float): Total simulation time.
        nx (int): Number of spatial steps.
        nt (int): Number of time steps.
        D (float): Diffusivity coefficient of the medium.
        velocity (float): Advection velocity.

    Returns:
        tuple: A tuple containing:
            - r_diff (float): Diffusion accuracy factor (D * Δt / Δx²).
            - r_adv (float): Advection accuracy factor (velocity * Δt / Δx).
    """
    dx, dt = calculate_discretization(L, T, nx, nt)
    r_diff = D * dt / dx**2
    r_adv = velocity * dt / dx

    return r_diff, r_adv

def check_accuracy_guidelines(L, T, nx, nt, D, velocity):
    """
    Validate the stability conditions for accuracy in the 1D advection-diffusion equation.

    Ensures that the diffusion stability factor (r_diff) is below 0.5 and the 
    advection stability factor (r_adv) is below 1.0 for numerical accuracy.
    This does not enforce strict stability limits but serves as a guideline for
    maintaining solution accuracy and physical interpretability.

    Args:
        L (float): Length of the domain.
        T (float): Total simulation time.
        nx (int): Number of spatial steps.
        nt (int): Number of time steps.
        D (float): Diffusivity coefficient of the medium.
        velocity (float): Advection velocity.

    Warns:
        UserWarning: If the accuracy conditions are not met.
    """
    r_diff, r_adv = calculate_accuracy_factors(L, T, nx, nt, D, velocity)
    
    if r_diff > 0.5:
        warnings.warn(f"Potential accuracy issue: r_diff={r_diff}. Recommended r_diff < 0.5 for accuracy.", UserWarning)
    if r_adv > 1.0:
        warnings.warn(f"Potential accuracy issue: r_adv={r_adv}. Recommended r_adv < 1.0 for accuracy.", UserWarning)
        
def setup_gaussian_pulse(L, nx, x0=None, sigma=None):
    """
    Generates an initial Gaussian pulse distribution on a 1D grid.

    Args:
        L (float): Length of the domain.
        nx (int): Number of spatial steps.
        x0 (float, optional): Center of the Gaussian pulse. Defaults to L/2.
        sigma (float, optional): Standard deviation (controls the width of the pulse). 
                                 Defaults to L/20.

    Returns:
        numpy.ndarray: Initial concentration profile.
    """
    if x0 is None:
        x0 = L / 2  
    if sigma is None:
        sigma = L / 20  
        
    dx, dt = calculate_discretization(L, T=1.0, nx=nx, nt=5) # Use dummy T and any nt > 2
    x = np.linspace(0, (nx - 1) * dx, nx)
    
    return np.exp(-0.5 * ((x - x0) / sigma)**2)

def create_matrices(nx, r_diff, r_conv):
    """
    Creates matrices A and B for the Crank-Nicolson method used in solving 
    the 1D advection-diffusion equation.
 
    Args:
        nx (int): Number of spatial grid points.
        r_diff (float): Diffusion number (D * dt / dx^2).
        r_conv (float): Convection number (v * dt / dx).

    Returns:
        tuple of numpy.ndarray: Matrices A and B used in the Crank-Nicolson 
        time-stepping scheme, where A is the left-hand side matrix and B is 
        the right-hand side matrix.
    """
    # Main diagonals
    A_main = np.ones(nx) * (1 + r_diff)
    B_main = np.ones(nx) * (1 - r_diff)

    # Upper and lower diagonals
    A_upper = np.ones(nx - 1) * (-0.5 * r_diff - 0.25 * r_conv)  
    A_lower = np.ones(nx - 1) * (-0.5 * r_diff + 0.25 * r_conv)  

    B_upper = np.ones(nx - 1) * (0.5 * r_diff + 0.25 * r_conv) 
    B_lower = np.ones(nx - 1) * (0.5 * r_diff - 0.25 * r_conv)  

    # Construct matrices
    A = np.diag(A_main) + np.diag(A_upper, k=1) + np.diag(A_lower, k=-1)
    B = np.diag(B_main) + np.diag(B_upper, k=1) + np.diag(B_lower, k=-1)

    return A, B