# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

#=
This example demonstrates linearized model accuracy by comparing:
1. Nonlinear SymbolicAWEModel model simulation 
2. Linearized state-space model simulation

Both models start from the same operating point and are subjected
to identical steering inputs. The resulting state trajectories are
plotted together to visualize how well the linearized model
approximates the nonlinear dynamics.
=#

using Timers
tic()
@info "Loading packages "

PLOT = true
if PLOT
    using Pkg
    if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, LaTeXStrings, ControlSystemsBase
end

using KiteModels, LinearAlgebra, Statistics, OrdinaryDiffEqCore 
using ModelingToolkit
using ModelingToolkit: t_nounits
toc()

# TODO: use sparse autodiff

# Simulation parameters
dt = 0.001
total_time = 1.0  # Increased from 0.1s to 1.0s for better dynamics observation
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 5.0      # Magnitude of steering input [Nm]

# Initialize model
set = load_settings("system_ram.yaml")
set.segments = 3
set_values = [-50.0, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "simple_ram"

@info "Creating SymbolicAWEModel model..."
s = SymbolicAWEModel(set)
s.set.abs_tol = 1e-2
s.set.rel_tol = 1e-2
toc()

# Define outputs for linearization - heading
lin_outputs = @variables heading(t_nounits)[1]

# Initialize at elevation with linearization outputs
s.sys_struct.winches[2].tether_length += 0.2
s.sys_struct.winches[3].tether_length += 0.2
KiteModels.init_sim!(s; 
    remake=false,
    reload=true,
    lin_outputs  # Specify which outputs to track in linear model
)
sys = s.sys

@show rad2deg(s.integrator[sys.elevation[1]])


@info "System initialized at:"
toc()

# --- Stabilize system at operating point ---
@info "Stabilizing system at operating point..."
s.integrator.ps[sys.stabilize] = true
stabilization_steps = Int(10 ÷ dt)
for i in 1:stabilization_steps
    next_step!(s; dt, vsm_interval=0.05÷dt)
end
s.integrator.ps[sys.stabilize] = false

# --- Linearize at operating point ---
@info "Linearizing system at operating point..."
@time (; A, B, C, D) = KiteModels.linearize(s)
@time (; A, B, C, D) = KiteModels.linearize(s)
@show norm(A)
@info "System linearized with matrix dimensions:" A=size(A) B=size(B) C=size(C) D=size(D)

sys = ss(A,B,C,D)
bode_plot(sys[1,1]; from=1e-4)
bode_plot(sys[1,2]; from=1e-4)
bode_plot(sys[1,3]; from=1e-4)
