using Revise, KiteModels, OrdinaryDiffEq
dt = 0.001
time = 0.01
s = KPS4_3L(KCU(se()))
# integrator = init_sim!(s; modeling_toolkit=true)
pos, vel = init_pos_vel(s)
simple_sys, sys = KiteModels.model!(s, pos, vel);
tspan = (0.0, dt)
prob = ODEProblem(simple_sys, nothing, tspan)
tol=1e-6
integrator = init(prob, TRBDF2(); dt=dt, abstol=tol, reltol=tol, save_everystep=false)
for (j, time) in pairs(0:dt:time)
    @time step!(integrator, dt, true)
end
nothing

# There are 1420 highest order derivative variables and 1374 equations