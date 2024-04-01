using Sundials, ControlPlots

G_EARTH  = [0.0, 0.0, -9.81] # gravitational acceleration
  
# Example one: Falling mass
# State vector y   = mass.pos, mass.vel
# Derivative   yd  = mass.vel, mass.acc
# Residual     res = (y.vel - yd.vel), (yd.acc - G_EARTH)     
function res!(res, yd, y, p, time)
    res[1:3] .= y[4:6] - yd[1:3]
    res[4:6] .= yd[4:6] - G_EARTH 
end

function run_example()
    # Set the initial conditions
    t_final = 10.0              # Final simulation time
    vel_0 = [0.0, 0.0, 50.0]    # Initial velocity
    pos_0 = [0.0, 0.0,  0.0]    # Initial position
    acc_0 = [0.0, 0.0, -9.81]   # Initial acceleration
    y0 = append!(pos_0, vel_0)  # Initial pos, vel
    yd0 = append!(vel_0, acc_0) # Initial vel, acc

    differential_vars = ones(Bool, length(y0))
    solver  = IDA(linear_solver=:GMRES, max_order = 4)
    tspan   = (0.0, t_final) 
    abstol  = 0.0006 # max error in m/s and m
    s = nothing

    prob    = DAEProblem(res!, yd0, y0, tspan, s, differential_vars=differential_vars)
    integrator = Sundials.init(prob, solver, abstol, reltol=0.001)

    dt = 0.05
    pos_z=Float64[]
    vel_z=Float64[]
    for t in 0:0.05:t_final    
        Sundials.step!(integrator, dt, true)
        push!(pos_z, integrator.u[3])
        push!(vel_z, integrator.u[6])
    end
    
    # plot the result

    # plt.ax1 = plt.subplot(111) 
    # plt.ax1.set_xlabel('time [s]')
    # plt.plot(time, pos_z, color="green")
    # plt.ax1.set_ylabel('pos_z [m]')  
    # plt.ax1.grid(True) 
    # plt.ax2 = plt.twinx()  
    # plt.ax2.set_ylabel('vel_z [m/s]')   
    # plt.plot(time, vel_z, color="red")    
    # plt.show()
    integrator
end
