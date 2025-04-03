# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:56 2025

@author: madhu
"""

import pytest
from math import isclose
import warnings
import numpy as np
from functions import (calculate_discretization, calculate_and_check_accuracy_factors,
                       setup_gaussian_pulse, create_matrices, apply_boundary_conditions, solve_advection_diffusion_CN,
                       solve_advection_diffusion_analytical)

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
    
    GIVEN: A set of parameters for domain length, total time, and number of spatial and temporal grid points.
    WHEN: The calculate_discretization function is called.
    THEN: The computed dx and dt should be close to L/(nx-1) and T/(nt-1) respectively.
    
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
    Test that calculate_discretization raises a ValueError when nx and / or nt is less than 2.
    
    GIVEN: Invalid nx and / or nt (less than 2).
    WHEN: The calculate_discretization function is called.
    THEN: A ValueError should be raised.
    
    """
    with pytest.raises(ValueError):
        calculate_discretization(L, T, nx, nt)

@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_calculate_and_check_accuracy_factors(params):
    """
    Test that calculate_and_check_accuracy_factors computes the correct accuracy factors.
    
    GIVEN: A set of parameters for domain length, total time, number of spatial grid points
           number of temporal grid points, diffusivity, and velocity.
    WHEN: The calculate_and_check_accuracy_factors function is called.
    THEN: The computed r_diff and r_adv should be close to the manually calculated expected values.
    
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
    r_diff, r_adv = calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)
    
    assert isclose(r_diff, r_diff_expected, rel_tol=1e-9), \
        f"Diffusion factor mismatch: expected {r_diff_expected}, got {r_diff}"
    assert isclose(r_adv, r_adv_expected, rel_tol=1e-9), \
        f"Advection factor mismatch: expected {r_adv_expected}, got {r_adv}"
        
@pytest.mark.parametrize("params", unstable_adv_diffusion_cases)
def test_calculate_and_check_accuracy_factors_warns(params):
    """
    Test that calculate_and_check_accuracy_factors raises a warning for unstable configurations.
    
    GIVEN: A set of numerically unstable parameters for the advection-diffusion equation.
    WHEN: check_accuracy_guidelines is called.
    THEN: Warning should be raised indicating instability and accuracy issues.
    
    """
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")  
        calculate_and_check_accuracy_factors(**params)

@pytest.mark.parametrize("L, T, nx, nt, D, velocity", [
    (1.0, 0.1, 1, 11, 0.01, 0.5),  # nx is less than 2
    (1.0, 0.1, 11, 1, 0.01, 0.5),  # nt is less than 2
    (1.0, 0.1, 0, 0, 0.01, 0.5),   # both nx and nt are invalid
])
def test_calculate_and_check_accuracy_factors_invalid_inputs(L, T, nx, nt, D, velocity):
    """
    Test that calculate_and_check_accuracy_factors raises ValueError when invalid nx/nt are passed.
    
    GIVEN: Invalid nx and / or nt (less than 2) along with other valid parameters.
    WHEN: calculate_and_check_accuracy_factors is called.
    THEN: A ValueError is propagated from calculate_discretization.
    
    """
    with pytest.raises(ValueError):
        calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)

@pytest.mark.parametrize("params", advection_diffusion_cases)
def test_default_setup_gaussian_pulse(params):
    """
    Test that the setup_gaussian_pulse generates a Gaussian pulse with correct properties for default x0 and sigma.
    
    GIVEN: A set of parameters describing the domain length (L), number of spatial steps (nx).
    WHEN: The setup_gaussian_pulse function is called with these parameters.
    THEN: The resulting Gaussian pulse should be a numpy array with the correct shape, and the values 
          of the pulse should be close to that of a gaussian pulse with default x0 and sigma.

    """
    L = params["L"]
    nx = params["nx"]
    x0_default = L / 3
    sigma_default = L / 10
    
    result = setup_gaussian_pulse(L, nx, x0 = x0_default, sigma = sigma_default)  
   
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
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, and velocity
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
    Test that the create_matrices generates Crank-Nicolson matrices A and B with correct properties.
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, and velocity.
    WHEN: The create_matrices function is called with computed diffusion and advection factors.
    THEN: The returned matrices A and B should have shape (nx, nx) and contain expected values along diagonals,
          as required for the Crank-Nicolson method for the 1-D advection-diffusion equation.
          
    """
    L, T, D, velocity, nx, nt = params["L"], params["T"], params["D"], params["velocity"], params["nx"], params["nt"] 
    r_diff, r_adv = calculate_and_check_accuracy_factors(L, T, nx, nt, D, velocity)
    A, B = create_matrices(nx, r_diff, r_adv)

    assert A.shape == (nx, nx), f"A shape mismatch: {A.shape}"
    assert B.shape == (nx, nx), f"B shape mismatch: {B.shape}"
    assert np.allclose(np.diag(A), 1 + r_diff), "A main diagonal incorrect"
    assert np.allclose(np.diag(B), 1 - r_diff), "B main diagonal incorrect"
    assert np.allclose(np.diag(A, k=1), -0.5 * r_diff + 0.25 * r_adv), "A upper diagonal incorrect"
    assert np.allclose(np.diag(A, k=-1), -0.5 * r_diff - 0.25 * r_adv), "A lower diagonal incorrect"
    assert np.allclose(np.diag(B, k=1), 0.5 * r_diff - 0.25 * r_adv), "B upper diagonal incorrect"
    assert np.allclose(np.diag(B, k=-1), 0.5 * r_diff + 0.25 * r_adv), "B lower diagonal incorrect"

@pytest.mark.parametrize("nx", [2, 3, 5])
@pytest.mark.parametrize("r_diff", [0.0, 0.1])
@pytest.mark.parametrize("r_adv", [0.0, 0.2])
def test_create_matrices_edge_cases(nx, r_diff, r_adv):
    """
    Test create_matrices across edge cases involving zero and non-zero r_diff and r_adv values with small nx values.

    GIVEN: Combinations of small nx values with zero and non-zero r_diff and r_adv values.
    WHEN: create_matrices is called with each combination.
    THEN: The returned matrices A and B should have shape (nx, nx) and contain expected values along diagnols,
          as required for the Crank-Nicolson method for the 1-D advection-diffusion equation.
          
    """
    A, B = create_matrices(nx, r_diff, r_adv)

    assert A.shape == (nx, nx), f"A shape mismatch: {A.shape}"
    assert B.shape == (nx, nx), f"B shape mismatch: {B.shape}"
    assert np.allclose(np.diag(A), 1 + r_diff), "A main diagonal incorrect"
    assert np.allclose(np.diag(B), 1 - r_diff), "B main diagonal incorrect"
    assert np.allclose(np.diag(A, k=1), -0.5 * r_diff + 0.25 * r_adv), "A upper diagonal incorrect"
    assert np.allclose(np.diag(A, k=-1), -0.5 * r_diff - 0.25 * r_adv), "A lower diagonal incorrect"
    assert np.allclose(np.diag(B, k=1), 0.5 * r_diff - 0.25 * r_adv), "B upper diagonal incorrect"
    assert np.allclose(np.diag(B, k=-1), 0.5 * r_diff + 0.25 * r_adv), "B lower diagonal incorrect"
    
@pytest.mark.parametrize("shape", [(3, 3), (5, 5), (10, 10)])
def test_boundary_rows_modified_correctly(shape):
    """
    GIVEN: A and B matrices of different sizes with all elements equal to one.
    WHEN: apply_boundary_conditions is called.
    THEN: First and last rows should be zeroed out, except (0,0) and (-1,-1) which should be 1.
    
    """
    assert len(shape) == 2, "Shape must be 2-dimensional"
    
    A = np.ones(shape)
    B = np.ones(shape)

    A_result, B_result = apply_boundary_conditions(A.copy(), B.copy())

    for M in [A_result, B_result]:
        assert np.all(M[0, 1:] == 0), "Top row not zeroed properly"
        assert M[0, 0] == 1, "Top-left corner not set to 1"
        assert np.all(M[-1, :-1] == 0), "Bottom row not zeroed properly"
        assert M[-1, -1] == 1, "Bottom-right corner not set to 1"

@pytest.mark.parametrize("shape", [(3, 3), (5, 5), (10, 10)])
@pytest.mark.parametrize("matrix_func", ["eye", "zeros", "rand"])
def test_inner_values_preserved(shape, matrix_func):
    """
    GIVEN: A and B matrices initialized with various inner values.
    WHEN: apply_boundary_conditions is called.
    THEN: The interior of the matrices should remain unchanged.
    
    """
    assert len(shape) == 2, "Shape must be 2-dimensional"
    
    def generate_matrix(func_name, shape):
        if func_name == "eye":
            return np.eye(N=shape[0], M=shape[1])
        elif func_name == "zeros":
            return np.zeros(shape)
        elif func_name == "rand":
            return np.random.rand(*shape)
        else:
            raise ValueError("Unknown matrix function")

    A = generate_matrix(matrix_func, shape)
    B = generate_matrix(matrix_func, shape)

    A_inner_expected = A[1:-1, 1:-1].copy()
    B_inner_expected = B[1:-1, 1:-1].copy()

    A_result, B_result = apply_boundary_conditions(A.copy(), B.copy())

    assert np.allclose(A_result[1:-1, 1:-1], A_inner_expected), "Interior of A changed unexpectedly"
    assert np.allclose(B_result[1:-1, 1:-1], B_inner_expected), "Interior of B changed unexpectedly"
    
@pytest.mark.parametrize(
    "params",
    [
        {**case, "x0": case["L"] / 3, "sigma": case["L"] / 10}
        for case in advection_diffusion_cases
    ]
)
def test_solve_advection_diffusion_CN_stability(params):
    """
    Test that solve_advection_diffusion_CN behaves normally with configurations satisfying accuracy guidelines.
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, and velocity.
    WHEN: solve_advection_diffusion_CN is executed.
    THEN: The solution should remain stable, boundary conditions should be applied and have smooth time solution.
    
    """
    nx = params['nx']
    nt = params['nt']
    x, u = solve_advection_diffusion_CN(**params)
    max_time_step_change = np.max(np.abs(u[:, 1:] - u[:, :-1]))
    assert np.all(np.isfinite(u)), "Solution contains NaN or Inf values."
    assert np.allclose(u[0, :], 0), "Left boundary condition not satisfied"
    assert np.allclose(u[-1, :], 0), "Right boundary condition not satisfied"
    assert max_time_step_change < 1.0, f"Instability detected: max change {max_time_step_change} too high."
    assert u.shape == (nx, nt), f"Expected shape {(nx, nt)}, got {u.shape}"
    

@pytest.mark.parametrize(
    "params",
    [
        {**case, "x0": case["L"] / 3, "sigma": case["L"] / 10}
        for case in unstable_adv_diffusion_cases
    ]
)
def test_solve_advection_diffusion_CN_warns(params):
    """
    Test that solve_advection_diffusion_CN raises a warning for configurations not satisfying accuracy guidelines.
    
    GIVEN: A set of numerically unstable parameters for the advection-diffusion equation.
    WHEN: solve_advection_diffusion_CN is executed.
    THEN: A warning should be raised indicating instability and accuracy issues.
    
    """
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        solve_advection_diffusion_CN(**params)

@pytest.mark.parametrize(
    "params",
    [
        {**param, "x0": param["L"] / 3, "sigma": param["L"] / 10, "num_reflections": 5} 
        for param in advection_diffusion_cases
    ],
)
def test_solve_advection_diffusion_analytical(params):
    """
    Test that solve_advection_diffusion_analytical behaves correctly with same number of reflections,
    initial conditions, and pulse widths.

    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, velocity and initial conditions.
    WHEN: solve_advection_diffusion_analytical is executed.
    THEN: The solution should remain stable and have a smooth time evolution.
    
    """
    nx = params['nx']
    nt = params['nt']
    x, u = solve_advection_diffusion_analytical(**params)
    max_change = np.max(np.abs(np.diff(u, axis=1)))
    assert np.all(np.isfinite(u)), "Solution contains NaN or Inf values."   
    assert max_change < 1.0, f"Instability detected: max change {max_change} too high."
    assert u.shape == (nx, nt), f"Expected shape {(nx, nt)}, got {u.shape}"

@pytest.mark.parametrize(
    "params, x0, sigma, num_reflections", 
    [(param, x0, sigma, num_reflections) for param in advection_diffusion_cases 
     for x0, sigma in [(5.0, 1.0), (3.0, 0.5), (7.0, 2.0)] 
     for num_reflections in [1, 5, 10, 20, -1, -5, -10]]
)
def test_solve_advection_diffusion_analytical_reflections(params, x0, sigma, num_reflections):
    """
    Test that solve_advection_diffusion_analytical behaves correctly with different numbers of reflections,
    initial conditions, and pulse widths.
    
    GIVEN: A set of parameters for domain length, total time, discretization sizes, diffusivity, velocity,
           initial conditions (x0, sigma), and number of reflections.
    WHEN: solve_advection_diffusion_analytical is called.
    THEN: If num_reflections < 0, a ValueError is raised; otherwise, the function executes, producing a stable and smooth time evolution solution.
    
    """
    if num_reflections < 0:
        with pytest.raises(ValueError):
            solve_advection_diffusion_analytical(num_reflections=num_reflections, x0=x0, sigma=sigma, **params)
    else:
        nx = params['nx']
        nt = params['nt']
        x, u = solve_advection_diffusion_analytical(num_reflections=num_reflections, x0=x0, sigma=sigma, **params)
        max_change = np.max(np.abs(np.diff(u, axis=1)))
        assert np.all(np.isfinite(u)), "Solution contains NaN or Inf values."   
        assert max_change < 1.0, f"Instability detected: max change {max_change} too high."
        assert u.shape == (nx, nt), f"Expected shape {(nx, nt)}, got {u.shape}"