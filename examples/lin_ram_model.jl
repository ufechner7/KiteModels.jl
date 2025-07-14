# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# This example demonstrates how to:
# - Load a system configuration from a YAML file,
# - Initialize a SymbolicAWEModel with specified winch torques,
# - Stabilize the system at its operating point,
# - Linearize the system at that point,
# - And plot the Bode plots of the resulting linearized system.

using Timers
tic()
@info "Loading packages "

PLOT = true
using Pkg
if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings, ControlSystemsBase

using SymbolicAWEModels, KiteUtils, LinearAlgebra, Statistics
using ModelingToolkit
using ModelingToolkit: t_nounits
toc()

# Initialize model
set = Settings("system_ram.yaml")
set_values = [-50.0, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]

@info "Creating SymbolicAWEModel..."
s = SymbolicAWEModel(set)
toc()

# Define outputs for linearization - heading
@variables begin
    heading(t_nounits)[1]
    angle_of_attack(t_nounits)[1]
    tether_len(t_nounits)[1:3]
    winch_force(t_nounits)[1:3]
end
lin_outputs = [heading[1], angle_of_attack[1], tether_len[1], winch_force[1]]
@info "Linear outputs: $lin_outputs"

# Initialize at elevation with linearization outputs
SymbolicAWEModels.init!(s; lin_outputs)
sys = s.sys

@info "System initialized at:"
toc()

# --- Stabilize system at operating point ---
@info "Stabilizing system at operating point..."
SymbolicAWEModels.find_steady_state!(s)

# --- Linearize at operating point ---
@info "Linearizing system at operating point..."
(; A, B, C, D) = SymbolicAWEModels.linearize!(s)
@time (; A, B, C, D) = SymbolicAWEModels.linearize!(s)
@info "System linearized with matrix dimensions:" A=size(A) B=size(B) C=size(C) D=size(D)

sys = ss(A,B,C,D)
bode_plot(sys[1,1]; from=1e-4)
bode_plot(sys[1,2]; from=1e-4)
bode_plot(sys[1,3]; from=1e-4)
