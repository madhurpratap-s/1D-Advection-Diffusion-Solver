# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:28 2025

@author: madhu
"""

import configparser
from functions import calculate_and_check_accuracy_factors

def process_configuration_file(config_file):
    """
    Reads and processes the configuration file for 1-D Advection-Diffusion solver.

    Args:
        config_file (str): The path to the configuration file in INI format. The file should contain the following sections:
            - [simulation_paramters]: Contains the simulation parameters.
                - length (float): Length of the spatial domain.
                - nx (int): Number of spatial discretization points.
                - total time (float): Total time for the simulation.
                - nt (int): Number of time steps for the simulation.
                - alpha (float): Diffusivity coefficient.
                - velocity (float): Advection velocity.
            - [paths]: Specifies file paths for storing solutions.
                - numerical_solution (str): Path to save the numerical solution (in .npy format).
                - analytical_solution (str): Path to save the analytical solution (in .npy format).

    Process:
        1. Reads the configuration file and extracts the required simulation parameters and output file paths.
        2. Calculates and checks if the accuracy factors r_diff and r_adv meet accuracy guidelines.
        3. Computes both numerical and analytical solutions based on the provided parameters.
        4. Saves the computed solutions to the specified paths.
        5. Displays plots comparing the numerical and analytical solutions.
    
    Warns:
        UserWarning: If the accuracy conditions are not met (if r_diff > 0.5 or r_adv > 1.0).

    Raises:
        ValueError: If nx or nt are less than 2, as they must represent valid discretization parameters.
        ValueError: If num_reflections is less than 1 in solve_advection_diffusion_analytical.
    """

    config = configparser.ConfigParser()
    config.read(config_file)

    # Read simulation parameters from the config file
    L = float(config.get('simulation_paramters', 'length'))
    T = float(config.get('simulation_paramters', 'total time'))
    nx = int(config.get('simulation_paramters', 'nx'))
    nt = int(config.get('simulation_paramters', 'nt'))
    D = float(config.get('simulation_paramters', 'diffusivity'))
    velocity = float(config.get('simulation_paramters', 'velocity'))
    
    calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)

    # Read output paths
    numerical_solution_path = config.get('paths', 'numerical_solution')
    analytical_solution_path = config.get('paths', 'analytical_solution')

    print(f"Running simulation with nx={nx}, nt={nt}")