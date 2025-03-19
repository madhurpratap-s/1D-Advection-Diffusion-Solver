# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:56 2025

@author: madhu
"""

import pytest
from math import isclose
import warnings
import numpy as np
from functions import calculate_discretization, calculate_accuracy_factors, check_accuracy_guidelines, setup_gaussian_pulse, create_matrices

# Numerical test cases for different configurations
advection_diffusion_cases = [
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 0.0, "nx": 21, "nt": 41},    # pure diffusion, low res
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 1.0, "nx": 21, "nt": 41},    # advection-diffusion
    {"L": 0.5, "T": 0.1, "D": 0.05, "velocity": 0.5, "nx": 51, "nt": 101},   # short rod, moderate flow, high resolution
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 0.0, "nx": 41, "nt": 81},     # diffusion only, medium domain
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 2.0, "nx": 101, "nt": 401},   # high advection, fine grid
    {"L": 3.0, "T": 0.5, "D": 1e-4, "velocity": 1.0, "nx": 151, "nt": 301},  # dominant advection, low diffusion
    {"L": 5.0, "T": 1.0, "D": 0.5, "velocity": 0.2, "nx": 101, "nt": 401},   # long rod, high diffusion
    {"L": 5.0, "T": 1.0, "D": 0.1, "velocity": 0.5, "nx": 51, "nt": 101},    # long rod, moderate everything
    {"L": 1.0, "T": 0.1, "D": 0.01, "velocity": 0.0, "nx": 11, "nt": 11},    # low resolution baseline
]

unstable_adv_diffusion_cases = [
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 2.0, "nx": 101, "nt": 201},  # high advection, fine grid
    {"L": 5.0, "T": 1.0, "D": 0.5, "velocity": 3.0, "nx": 51, "nt": 101},   # both high diffusion and advection
    {"L": 1.0, "T": 0.1, "D": 1.0, "velocity": 15.0, "nx": 11, "nt": 11},   # very high diffusion and advection
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
        
@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_check_accuracy_guidelines_no_warning(params):
    """
    Test that check_accuracy_guidelines does not raise a warning for stable configurations for accuracy.
    
    GIVEN: A set of numerically stable parameters for the advection-diffusion equation.
    WHEN: check_accuracy_guidelines is called with these parameters.
    THEN: No warnings should be raised, indicating stability and accuracy conditions are met.
    
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")  
        check_accuracy_guidelines(**params)
        
@pytest.mark.parametrize("params", unstable_adv_diffusion_cases)
def test_check_accuracy_guidelines_warns(params):
    """
    Test that check_accuracy_guidelines raises a warning for unstable configurations.
    
    GIVEN: A set of numerically unstable parameters for the advection-diffusion equation.
    WHEN: check_accuracy_guidelines is called with these parameters.
    THEN: Warnings should be raised indicating instability and accuracy issues.
    
    """
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")  
        check_accuracy_guidelines(**params)

@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_default_setup_gaussian_pulse(params):
    """
    Test that the setup_gaussian_pulse generates a Gaussian pulse with correct propertie for default x0 and sigma.
    
    GIVEN: A set of parameters describing the domain length (L), number of spatial steps (nx).
    WHEN: The setup_gaussian_pulse function is called with these parameters.
    THEN: The resulting Gaussian pulse should be a numpy array with the correct shape, and the values 
          of the pulse should be close to that of a gaussian pulse with default x0 and sigma.

    """
    L = params["L"]
    nx = params["nx"]
    result = setup_gaussian_pulse(L, nx)  
   
    x0_default = L / 2
    sigma_default = L / 20
   
    x = np.linspace(0, L, nx)
    expected = np.exp(-0.5 * ((x - x0_default) / sigma_default)**2)
   
    assert isinstance(result, np.ndarray), "Result should be a numpy array"
    assert result.shape == (nx,), f"Expected shape {(nx,)}, got {result.shape}"
    np.testing.assert_allclose(result, expected, atol=1e-7,
             err_msg="Gaussian pulse does not match expected values with default x0 and sigma")
    
@pytest.mark.parametrize("params", advection_diffusion_cases)
@pytest.mark.parametrize("x0, sigma", [(0.3, 0.05), (0.7, 0.1), (1.0, 0.2)])
def test_custom_setup_gaussian_pulse(params, x0, sigma):
    """
    Test that the setup_gaussian_pulse generates a Gaussian pulse with correct properties for custom x0 and sigma.
    
    GIVEN: A set of parameters describing the domain length (L), number of spatial steps (nx),
           along with custom values for the center of the pulse (x0) and the width (sigma).
    WHEN: The setup_gaussian_pulse function is called with these parameters.
    THEN: The resulting Gaussian pulse should be a numpy array with the correct shape, and the values
          of the pulse should be close to that of a gaussian pulse with custom values of x0 and sigma.

    """
    L = params["L"]
    nx = params["nx"]
    result = setup_gaussian_pulse(L, nx, x0=x0, sigma=sigma)  
    
    x = np.linspace(0, L, nx)
    expected = np.exp(-0.5 * ((x - x0) / sigma)**2)
    
    assert isinstance(result, np.ndarray), "Result should be a numpy array"
    assert result.shape == (nx,), f"Expected shape {(nx,)}, got {result.shape}"
    np.testing.assert_allclose(result, expected, atol=1e-7, 
             err_msg="Gaussian pulse does not match expected values for given x0 and sigma")

@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_create_matrices(params):
    """
    Test that the create_matrices generates Crank-Nicolson matrices A and B with correct structure.
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, and velocity.
    WHEN: The create_matrices function is called with computed diffusion and advection factors.
    THEN: The returned matrices A and B should have shape (nx, nx) and contain expected values along diagnols,
          as required for the Crank-Nicolson method for the 1-D advection-diffusion equation.
          
    """
    L, T, D, velocity, nx, nt = params["L"], params["T"], params["D"], params["velocity"], params["nx"], params["nt"] 
    r_diff, r_conv = calculate_accuracy_factors(L, T, nx, nt, D, velocity)
    A, B = create_matrices(nx, r_diff, r_conv)

    assert A.shape == (nx, nx), f"A shape mismatch: {A.shape}"
    assert B.shape == (nx, nx), f"B shape mismatch: {B.shape}"
    assert np.allclose(np.diag(A), 1 + r_diff), "A main diagonal incorrect"
    assert np.allclose(np.diag(B), 1 - r_diff), "B main diagonal incorrect"
    assert np.allclose(np.diag(A, k=1), -0.5 * r_diff - 0.25 * r_conv), "A upper diagonal incorrect"
    assert np.allclose(np.diag(A, k=-1), -0.5 * r_diff + 0.25 * r_conv), "A lower diagonal incorrect"
    assert np.allclose(np.diag(B, k=1), 0.5 * r_diff + 0.25 * r_conv), "B upper diagonal incorrect"
    assert np.allclose(np.diag(B, k=-1), 0.5 * r_diff - 0.25 * r_conv), "B lower diagonal incorrect"

@pytest.mark.parametrize("nx", [2, 3, 5])
@pytest.mark.parametrize("r_diff", [0.0, 0.1])
@pytest.mark.parametrize("r_conv", [0.0, 0.2])
def test_create_matrices_edge_cases(nx, r_diff, r_conv):
    """
    Test create_matrices across edge cases involving zero and non-zero r_diff and r_conv values with small nx values.

    GIVEN: Combinations of small nx values with zero and non-zero r_diff and r_conv values.
    WHEN: create_matrices is called with each combination.
    THEN: The returned matrices A and B should have shape (nx, nx) and contain expected values along diagnols,
          as required for the Crank-Nicolson method for the 1-D advection-diffusion equation.
    """
    A, B = create_matrices(nx, r_diff, r_conv)

    assert A.shape == (nx, nx), f"A shape mismatch: {A.shape}"
    assert B.shape == (nx, nx), f"B shape mismatch: {B.shape}"
    assert np.allclose(np.diag(A), 1 + r_diff), "A main diagonal incorrect"
    assert np.allclose(np.diag(B), 1 - r_diff), "B main diagonal incorrect"
    assert np.allclose(np.diag(A, k=1), -0.5 * r_diff - 0.25 * r_conv), "A upper diagonal incorrect"
    assert np.allclose(np.diag(A, k=-1), -0.5 * r_diff + 0.25 * r_conv), "A lower diagonal incorrect"
    assert np.allclose(np.diag(B, k=1), 0.5 * r_diff + 0.25 * r_conv), "B upper diagonal incorrect"
    assert np.allclose(np.diag(B, k=-1), 0.5 * r_diff - 0.25 * r_conv), "B lower diagonal incorrect"