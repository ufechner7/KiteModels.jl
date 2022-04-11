const SHOW_FRONT = false

function init2()
    kps4.set.alpha =  0.08163
    KiteModels.clear(kps4)
    kps4.set.l_tether = 150.0 # - kps4.set.height_k - kps4.set.h_bridle
    kps4.set.area = 10.18
    kps4.set.rel_side_area = 30.6
    kps4.set.mass = 6.21
    kps4.set.c_s = 0.6
    kps4.set.damping = 473.0     # unit damping coefficient
    kps4.set.c_spring = 614600.0 # unit spring coefficent
    kps4.set.width = 4.9622
    kps4.set.elevation = 70.7 
    kps4.set.profile_law = Int(EXPLOG)
    pos, vel = KiteModels.init_pos_vel(kps4)
    posd = copy(vel)
    veld = zero(vel)
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)
    kps4.stiffness_factor = 0.04
    KiteModels.init_springs(kps4)
    return pos, vel, posd, veld
end

function plot2d(x, z; zoom=1)
    xlabel = "x [m]"
    if SHOW_FRONT xlabel = "y [m]" end
    if zoom ==1
        x_max=maximum(x)
        z_max=maximum(z)
        plot(x,z, xlabel=xlabel, ylabel="z [m]", legend=false, xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5))
    else
        plot(x,z, xlabel=xlabel, ylabel="z [m]", legend=false)
    end
    plot!([x[7],x[10]],[z[7],z[10]], legend=false) #s6
    plot!([x[8],x[11]],[z[8],z[11]], legend=false) #s8
    plot!([x[9],x[11]],[z[9],z[11]], legend=false) #s7
    plot!([x[8],x[10]],[z[8],z[10]], legend=false) #s2
    plot!([x[7],x[11]] ,[z[7],z[11]],legend=false) #s5    
    plot!(x, z, seriestype = :scatter)
end

KiteModels.set_depower_steering(kps4, 0.25, 0.0)

x0 = Float64[] 
z0 = Float64[]
if typeof(kps4) <: KPS4
    pos, vel, posd, veld = init2()
    kps4.alpha_depower = deg2rad(2.2095658807330962) # from one point simulation
    kps4.set.alpha_zero = 0.0   
end

y0, yd0 = KiteModels.init(kps4)

for i in 1:length(kps4.pos)
     push!(x0, kps4.pos[i][1])
     push!(z0, kps4.pos[i][3])
end

find_steady_state(kps4, true)

x = Float64[] 
z = Float64[]

for i in 1:length(kps4.pos)
    if SHOW_FRONT
        push!(x, kps4.pos[i][2])
    else
        push!(x, kps4.pos[i][1])
    end
    push!(z, kps4.pos[i][3])
end

println("kite distance: $(norm(kps4.pos[end]))")
println(KiteModels.spring_forces(kps4))
println("alpha_depower [deg]: $(rad2deg(kps4.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")

integrator=KiteModels.init_sim(kps4, 1.0)

x1 = Float64[] 
z1 = Float64[]

for i in 1:length(kps4.pos)
    if SHOW_FRONT
        push!(x1, kps4.pos[i][2])
    else
        push!(x1, kps4.pos[i][1])
    end
    push!(z1, kps4.pos[i][3])
end

println("kite distance: $(norm(kps4.pos[end]))")
println(KiteModels.spring_forces(kps4))
println("alpha_depower [deg]: $(rad2deg(kps4.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")

#plot2d(x, z; zoom=1)
plot2d(x1, z1; zoom=1)

