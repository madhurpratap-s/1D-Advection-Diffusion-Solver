# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:52:29 2025

@author: madhu
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_1d_solutions(x, u_numerical, u_analytical, nt, T, L, nx, D, velocity):
    """
    Plots the initial, mid-time, and final solutions for both numerical and analytical 
    solutions of the 1D advection-diffusion equation.
    
    Args:
        x (numpy.ndarray): 1D array of spatial grid points.
        u_numerical (numpy.ndarray): 2D array (shape: [nx, nt]) containing the numerical solution.
        u_analytical (numpy.ndarray): 2D array (shape: [nx, nt]) containing the analytical solution.
        nt (int): Number of time steps.
        T (float): Total simulation time.
        L (float): Length of the domain.
        nx (int): Number of spatial discretization points.
        D (float): Diffusivity coefficient.
        velocity (float): Advection velocity.
        
    """
    indices = [0, nt // 2, nt - 1]  
    time_labels = ["t = 0", f"t = {T/2:.2f}", f"t = {T:.2f}"]  
    
    plt.figure(figsize=(10, 6))
    for idx, label in zip(indices, time_labels):
        plt.plot(x, u_numerical[:, idx], 'o-', label=f'Numerical {label}')
        plt.plot(x, u_analytical[:, idx], 'x--', label=f'Analytical {label}')
    
    plt.xlabel("Spatial coordinate, x")
    plt.ylabel("Solution, u")
    plt.title(f"1D Advection-Diffusion Solutions\n(Diffusivity = {D}, Velocity = {velocity})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_3d_surface(x, u_numerical, nt, T, L, nx, D, velocity):
    """
    Creates a 3D surface plot of the numerical solution over space and time.
    
    Args:
        x (numpy.ndarray): 1D array of spatial grid points.
        u_numerical (numpy.ndarray): 2D array (shape: [nx, nt]) containing the numerical solution.
        nt (int): Number of time steps.
        T (float): Total simulation time.
        L (float): Length of the domain.
        nx (int): Number of spatial discretization points.
        D (float): Diffusivity coefficient.
        velocity (float): Advection velocity.
        
    """
    # Create a time array corresponding to each time step
    t = np.linspace(0, T, nt)
    
    # Create a meshgrid for plotting
    T, X = np.meshgrid(t, x)
    
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, u_numerical, cmap='viridis', alpha=0.8)
    
    ax.set_xlabel('Spatial coordinate, x')
    ax.set_ylabel('Time, t')
    ax.set_zlabel('Solution, u')
    ax.set_title(f"Surface Plot of Numerical Solution\n(Diffusivity = {D}, Velocity={velocity})")
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.tight_layout()
    plt.show()
    