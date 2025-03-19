# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:52:20 2025

@author: madhu
"""

import warnings

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