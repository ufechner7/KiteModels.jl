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
    kps4.stiffness_factor = 0.5
    KiteModels.init_springs(kps4)
    return pos, vel, posd, veld
end

function plot2d(x, z; zoom=1)
    if zoom ==1
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false, xlims = (35.0, 55), ylims = (135.0, 155))
        plot!([x[7],x[10]],[z[7],z[10]], legend=false)
        plot!([x[8],x[11]],[z[8],z[11]], legend=false)
    else
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    end
    plot!(x, z, seriestype = :scatter)
end

KiteModels.set_depower_steering(kps4, 0.25, 0.0)

x0 = Float64[] 
z0 = Float64[]
if typeof(kps4) <: KPS4
    pos, vel, posd, veld = init2()
    kps4.alpha_depower = deg2rad(2.2095658807330962) # from one point simulation
    height = 134.14733504839947                      # from one point simulation
    kps4.stiffness_factor = 1.0
    kps4.set.alpha_zero = 0.0   
end

println(kps4.set.l_tether)
y0, yd0 = KiteModels.init(kps4)

for i in 1:length(kps4.pos)
     push!(x0, kps4.pos[i][1])
     push!(z0, kps4.pos[i][3])
end

find_steady_state(kps4, true)

x = Float64[] 
z = Float64[]

for i in 1:length(kps4.pos)
     push!(x, kps4.pos[i][1])
     push!(z, kps4.pos[i][3])
end

println("kite distance: $(norm(kps4.pos[end]))")
println(KiteModels.spring_forces(kps4))
println("alpha_depower [deg]: $(rad2deg(kps4.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")

plot2d(x, z; zoom=1)

