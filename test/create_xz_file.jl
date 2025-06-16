# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

using Timers
tic()
@info "Loading packages "

using KiteModels, LinearAlgebra, Statistics

toc()

# Simulation parameters
dt = 0.05
total_time = 10  # Longer simulation to see oscillations
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 10.0      # Magnitude of steering input [Nm]

# Initialize model
set = load_settings("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "ram"

@info "Creating wing, aero, vsm_solver, sys_struct and s:"
s = SymbolicAWEModel(set)
s.set.abs_tol = 1e-2
s.set.rel_tol = 1e-2
toc()

# Initialize at elevation
s.sys_struct.winches[2].tether_length += 0.2
s.sys_struct.winches[3].tether_length += 0.2
KiteModels.init_sim!(s; remake=false, reload=true)
sys = s.sys

@info "System initialized at:"
toc()
