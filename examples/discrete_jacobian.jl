#=
This example analyzes the input-output behavior of a ram air kite model by:

1. Creating a ram air kite model and initializing it at 60° elevation
2. Stabilizing the system
3. Testing the system response to different control inputs:
   - Left steering line torque
   - Right steering line torque

The script plots relationships between steering inputs and resulting angular velocities,
providing insight into the kite's steering behavior and control characteristics.
=#

using Timers
tic()
@info "Loading packages "

using KiteModels, LinearAlgebra, Statistics, OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqBDF
using DifferentiationInterface, FiniteDiff, SparseArrays, ADTypes
using ModelingToolkit: setu, getu, setp, getp
using ModelingToolkit

NOISE = 0.1
PLOT = true
if PLOT
    using Pkg
    if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
    using ControlPlots, LaTeXStrings
    import ControlPlots: plot
end
toc()


include(joinpath(@__DIR__, "plotting.jl"))

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
set.segments = 2
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = true
set.physical_model = "ram"

@info "Creating wing, aero, vsm_solver, point_system and s:"
s = RamAirKite(set)
s.set.abs_tol = 1e-3
s.set.rel_tol = 1e-3
toc()

# init_Q_b_w, R_b_w = KiteModels.measure_to_q(measure)
# init_kite_pos = init!(s.point_system, s.set, R_b_w)
# plot(s.point_system, 0.0; zoom=false, front=true)

measure = Measurement()

# Initialize at elevation
measure.sphere_pos .= deg2rad.([70.0 70.0; 1.0 -1.0])
KiteModels.init_sim!(s, measure; 
    adaptive=false, remake=false, reload=true, 
    solver=FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.8, max_iter=1000))
)
OrdinaryDiffEqCore.set_proposed_dt!(s.integrator, dt)
sys = s.sys

@info "System initialized at:"
toc()

sys = s.sys
# # Stabilize system
s.integrator.ps[sys.stabilize] = true
for i in 1:10÷dt
    next_step!(s; dt, vsm_interval=1)
end
s.integrator.ps[sys.stabilize] = false
# plot(s, 0.0; zoom=true, front=false)

# Function to step simulation with input u
function f(x, u, _, p)
    (s, set_x, set_u, get_x, dt) = p
    set_x(s.integrator, x)
    set_u(s.integrator, u)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=false)
    OrdinaryDiffEqCore.step!(s.integrator, dt)
    xnext = get_x(s.integrator)
    @show norm(xnext)
    return xnext
end

# Get initial state
x_vec = KiteModels.get_unknowns(s)
set_ix = setu(s.integrator, Initial.(x_vec)) # set_ix might not be needed anymore if prob is removed
set_x = setu(s.integrator, x_vec)
set_u = setu(s.integrator, collect(sys.set_values))
get_u = getu(s.integrator, collect(sys.set_values))
get_x = getu(s.integrator, x_vec)
x0 = get_x(s.integrator)
u0 = get_u(s.integrator)
p = (s, set_x, set_u, get_x, dt)

@show x0
x0 .+= (rand(length(x0)) .- 0.5) * NOISE
set_ix(s.integrator, x0)
OrdinaryDiffEqCore.reinit!(s.integrator)

f_x(x) = f(x, u0, nothing, p)
f_u(u) = f(x0, u, nothing, p)

x_idxs = Dict{Num, Int}()
for (idx, sym) in enumerate(x_vec)
    x_idxs[sym] = idx
end

m, n = length(x_vec), length(x_vec)
active_x = [
    [x_idxs[sys.free_twist_angle[i]] for i in 1:4]
    [x_idxs[sys.ω_b[i]] for i in 1:3]
    [x_idxs[sys.kite_pos[i]] for i in 1:3]
    [x_idxs[sys.kite_vel[i]] for i in 1:3]
]
# Collect all (row, col) pairs
rows = Int[]
cols = Int[]
for i in active_x, j in active_x
    push!(rows, i)
    push!(cols, j)
end
# Build the sparse boolean matrix
S = sparse(rows, cols, trues(length(rows)), m, n)


backend = AutoFiniteDiff()
sparse_backend = AutoSparse(
    backend;
    sparsity_detector=ADTypes.KnownJacobianSparsityDetector(S)
)
jac_prep = prepare_jacobian(f_x, backend, x0)
@time J = jacobian(f_x, jac_prep, backend, x0)
@time J = jacobian(f_x, jac_prep, backend, x0)
@show norm(J)
display(J)

