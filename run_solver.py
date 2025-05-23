# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:28 2025

@author: madhu
"""

import configparser
import argparse
import numpy as np
from functions import calculate_and_check_accuracy_factors, solve_advection_diffusion_CN, solve_advection_diffusion_analytical
from plot import plot_1d_solutions, plot_3d_surface

def process_configuration_file(config_file):
    """
    Reads and processes the configuration file for 1-D Advection-Diffusion solver.

    Args:
        config_file (str): The path or name of the configuration file in INI format. 
        If the file is in the same directory as the Python script, the file name is sufficient. 
        The file should contain the following sections:
            
            - [simulation_paramters]: Contains the simulation parameters.
                - length (float): Length of the spatial domain.
                - nx (int): Number of spatial discretization points.
                - total time (float): Total time for the simulation.
                - nt (int): Number of time steps for the simulation.
                - D (float): Diffusivity coefficient.
                - velocity (float): Advection velocity.
                - x0 (float, optional): Initial center of the Gaussian pulse. Defaults to length / 3.
                - sigma (float, optional): Standard deviation of the initial Gaussian pulse. Defaults to length / 10.
                - num_reflections (int, optional): Number of mirrored sources to consider. Defaults to 5.
                
            - [paths]: Specifies file paths for storing solutions.
                - numerical_solution (str): Path to save the numerical solution (in .npy format).
                - analytical_solution (str): Path to save the analytical solution (in .npy format).

    Process:
        1. Reads the configuration file and extracts the required simulation parameters and output file paths.
        2. Calculates and checks if the accuracy factors r_diff and r_adv meet accuracy guidelines.
        3. Computes both numerical and analytical solutions based on the provided parameters.
        4. Saves the computed solutions to the specified paths.
        5. Plots the 1D solutions of numerical and analytical methods with time together for comparison.
        6. Creates a 3D surface plot of the numerical solution over space and time.   
        
    Warns:
        UserWarning: If the accuracy conditions are not met (if r_diff > 0.5 or r_adv > 1.0).

    Raises:
        ValueError: If nx or nt are less than 2, as they must represent valid discretization parameters.
        ValueError: If num_reflections is less than 1 in solve_advection_diffusion_analytical.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
        
    L = float(config.get('simulation_parameters', 'length'))
    T = float(config.get('simulation_parameters', 'total time'))
    nx = int(config.get('simulation_parameters', 'nx'))
    nt = int(config.get('simulation_parameters', 'nt'))
    D = float(config.get('simulation_parameters', 'diffusivity'))
    velocity = float(config.get('simulation_parameters', 'velocity'))
    if 'num_reflections' in config['simulation_parameters']:
        num_reflections = int(config['simulation_parameters']['num_reflections'])
    else:
        num_reflections = 5  
    if 'x0' in config['simulation_parameters']:
        x0 = float(config['simulation_parameters']['x0'])
    else:
        x0 = None
    if 'sigma' in config['simulation_parameters']:
       sigma = float(config['simulation_parameters']['sigma'])
    else:
       sigma = None
    
    calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)

    numerical_solution_path = config.get('paths', 'numerical_solution')
    analytical_solution_path = config.get('paths', 'analytical_solution')

    print(f"Running simulation with nx={nx}, nt={nt}")
    
    x, u_numerical = solve_advection_diffusion_CN(L, T, nx, nt, D, velocity, x0, sigma)
    x, u_analytical = solve_advection_diffusion_analytical(L, T, nx, nt, D, velocity, x0, sigma, num_reflections)

    np.save(numerical_solution_path, u_numerical)
    np.save(analytical_solution_path, u_analytical)
    
    plot_1d_solutions(x, u_numerical, u_analytical, nt, T, L, nx, D, velocity)
    plot_3d_surface(x, u_numerical, nt, T, L, nx, D, velocity)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run 1-D advection-diffusion equation simulation using a configuration file."
    )
    parser.add_argument("config_file", nargs="?", default="configuration.txt")
    # The user can choose a specific configuration via command line or use the default one.
    args = parser.parse_args()
    process_configuration_file(args.config_file)
