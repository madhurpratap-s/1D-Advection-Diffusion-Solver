# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 21:53:56 2025

@author: madhu
"""

import pytest

# Numerical test cases for different configurations
advection_diffusion_cases = [
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 0.0, "nx": 21, "nt": 41},    # pure diffusion, low res
    {"L": 1.0, "T": 0.2, "D": 0.01, "velocity": 1.0, "nx": 21, "nt": 41},    # advection-diffusion
    {"L": 0.5, "T": 0.1, "D": 0.05, "velocity": 0.5, "nx": 51, "nt": 101},   # short rod, moderate flow, high resolution
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 0.0, "nx": 41, "nt": 81},     # diffusion only, medium domain
    {"L": 2.0, "T": 0.5, "D": 0.1, "velocity": 2.0, "nx": 101, "nt": 201},   # high advection, fine grid
    {"L": 3.0, "T": 0.5, "D": 1e-4, "velocity": 1.0, "nx": 151, "nt": 301},  # dominant advection, low diffusion
    {"L": 5.0, "T": 1.0, "D": 0.5, "velocity": 0.2, "nx": 101, "nt": 201},   # long rod, high diffusion
    {"L": 5.0, "T": 1.0, "D": 0.1, "velocity": 0.5, "nx": 51, "nt": 101},    # long rod, moderate everything
    {"L": 1.0, "T": 0.1, "D": 0.01, "velocity": 0.0, "nx": 11, "nt": 11},    # low resolution baseline
]
