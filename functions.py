# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:52:20 2025

@author: madhu
"""

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

