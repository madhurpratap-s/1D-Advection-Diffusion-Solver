# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:56 2025

@author: madhu
"""

import pytest
from math import isclose
import warnings
from functions import calculate_discretization, calculate_accuracy_factors, check_accuracy_guidelines

# Numerical test cases for different configurations
advection_diffusion_cases = [
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 0.0, "nx": 21, "nt": 41},    # pure diffusion, low res
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 1.0, "nx": 21, "nt": 41},    # advection-diffusion
    {"L": 0.5, "T": 0.1, "D": 0.05, "velocity": 0.5, "nx": 51, "nt": 101},   # short rod, moderate flow, high resolution
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 0.0, "nx": 41, "nt": 81},     # diffusion only, medium domain
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 2.0, "nx": 101, "nt": 401},   # high advection, fine grid
    {"L": 3.0, "T": 0.5, "D": 1e-4, "velocity": 1.0, "nx": 151, "nt": 301},  # dominant advection, low diffusion
    {"L": 5.0, "T": 1.0, "D": 0.5, "velocity": 0.2, "nx": 101, "nt": 401},    # long rod, high diffusion
    {"L": 5.0, "T": 1.0, "D": 0.1, "velocity": 0.5, "nx": 51, "nt": 101},    # long rod, moderate everything
    {"L": 1.0, "T": 0.1, "D": 0.01, "velocity": 0.0, "nx": 11, "nt": 11},    # low resolution baseline
]

@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_calculate_discretization(params):
    """
    Test that calculate_discretization computes the correct spatial and temporal step sizes.
    
    GIVEN: A set of parameters for domain length, total time, and discretization sizes.
    WHEN: The calculate_discretization function is called.
    THEN: The computed dx and dt should equal L/(nx-1) and T/(nt-1) respectively.
    
    """
    L = params["L"]
    T = params["T"]
    nx = params["nx"]
    nt = params["nt"] 
    
    # Calculate expected spatial and temporal steps
    expected_dx = L / (nx - 1)
    expected_dt = T / (nt - 1)

    dx, dt = calculate_discretization(L, T, nx, nt)
    assert isclose(dx, expected_dx, rel_tol=1e-9), f"Expected dx {expected_dx}, got {dx}"
    assert isclose(dt, expected_dt, rel_tol=1e-9), f"Expected dt {expected_dt}, got {dt}"

@pytest.mark.parametrize("L, T, nx, nt", [
    (1.0, 0.1, 1, 11),  # nx is less than 2
    (1.0, 0.1, 11, 1),  # nt is less than 2
    (1.0, 0.1, 0, 0)    # both nx or nt are invalid
])
def test_calculate_discretization_invalid_inputs(L, T, nx, nt):
    """
    Test that calculate_discretization raises a ValueError when nx or nt is less than 2.
    
    GIVEN: Invalid nx or nt (less than 2).
    WHEN: The calculate_discretization function is called.
    THEN: A ValueError should be raised.
    
    """
    with pytest.raises(ValueError):
        calculate_discretization(L, T, nx, nt)
        
@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_calculate_accuracy_factors(params):
    """
    Test that calculate_accuracy_factors computes the correct accuracy factors.
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, and velocity.
    WHEN: The calculate_accuracy_factors function is called.
    THEN: The computed r_diff and r_adv should match the manually calculated expected values.
    
    """
    L = params["L"]
    T = params["T"]
    nx = params["nx"]
    nt = params["nt"]
    D = params["D"]
    velocity = params["velocity"]
    
    dx_expected = L / (nx - 1)
    dt_expected = T / (nt - 1)
    
    # Compute expected accuracy factors
    r_diff_expected = D * dt_expected / dx_expected**2
    r_adv_expected = velocity * dt_expected / dx_expected
    
    # Get the actual factors from the function
    r_diff, r_adv = calculate_accuracy_factors(L, T, nx, nt, D, velocity)
    
    assert isclose(r_diff, r_diff_expected, rel_tol=1e-9), \
        f"Diffusion factor mismatch: expected {r_diff_expected}, got {r_diff}"
    assert isclose(r_adv, r_adv_expected, rel_tol=1e-9), \
        f"Advection factor mismatch: expected {r_adv_expected}, got {r_adv}"

@pytest.mark.parametrize("L, T, nx, nt, D, velocity", [
    (1.0, 0.1, 1, 11, 0.01, 0.5), # nx is less than 2
    (1.0, 0.1, 11, 1, 0.01, 0.5),  # nt is less than 2
    (1.0, 0.1, 0, 0, 0.01, 0.5), # both nx and nt are invalid
])
def test_calculate_accuracy_factors_invalid_inputs(L, T, nx, nt, D, velocity):
    """
    Test that calculate_accuracy_factors raises ValueError when invalid nx/nt are passed.
    
    GIVEN: Invalid nx or nt (less than 2) along with other valid parameters.
    WHEN: calculate_accuracy_factors is called.
    THEN: A ValueError is propagated from calculate_discretization.
    
    """
    with pytest.raises(ValueError):
        calculate_accuracy_factors(L, T, nx, nt, D, velocity)