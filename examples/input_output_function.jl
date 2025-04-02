using Revise, KiteModels, LinearAlgebra, VortexStepMethod
using ModelingToolkit
using ModelingToolkit: setu, getu

PLOT = false
if PLOT
    using ControlPlots
end

dt = 0.01
total_time = 1.5
vsm_interval = 5
steps = Int(round(total_time / dt))

set = se("system_3l.yaml")
set.segments = 2
set_values = [-50, -1.1, -1.1]

if !@isdefined s
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    s = RamAirKite(set, wing, aero, vsm_solver)
end
s.set.abs_tol = 1e-5
s.set.rel_tol = 1e-3
s.measure.sphere_pos .= deg2rad.([50.0 50.0; 1.0 -1.0])

# KiteModels.init_sim!(s)
@time KiteModels.reinit!(s)
sys = s.sys

sym_vec = KiteModels.get_unknowns_vec(s, s.point_system; non_observed=false)
set_x = setu(s.integrator, sym_vec)
get_x = getu(s.integrator, sym_vec)

# x is the non-stiff state, u is the inputs
function f!(x_plus, x, u, _, _)
    @time set_x(s.integrator, x)
    @time next_step!(s)
    x_plus .= get_x(s.integrator)
    nothing
end

x = get_x(s.integrator)
x_plus = copy(x)
@show x x_plus
f!(x_plus, x, zeros(3), nothing, nothing)
@show x x_plus
f!(x_plus, x, zeros(3), nothing, nothing)
@show x x_plus

nothing