function init_150()
    KiteModels.clear(kps4)
    kps4.set.l_tether = 150.0
    kps4.set.elevation = 70.0
    kps4.set.area = 10.18
    kps4.set.rel_side_area = 30.6
    kps4.set.v_wind = 9.1
    kps4.set.mass = 6.21
    kps4.set.c_s = 0.6
    kps4.set.damping = 473.0     # unit damping coefficient
    kps4.set.c_spring = 614600.0 # unit spring coefficent
    kps4.set.width = 4.9622
end

function init2()
    kps4.set.alpha = 1.0/7.0
    init_150()
    kps4.set.elevation = 70.7 
    kps4.set.profile_law = Int(EXP)
    kps4.set.l_tether = 150.0
    pos, vel = KiteModels.init_pos_vel(kps4)
    posd = copy(vel)
    veld = zero(vel)
    kps4.v_wind_gnd .= [9.1, 0.0, 0.0]
    height = 134.14733504839947
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)
    kps4.stiffness_factor = 0.5
    kps4.set.alpha = 1.0/7.0
    kps4.segment_length = kps4.set.l_tether/se().segments
    for i in 1:se().segments + KiteModels.KITE_PARTICLES + 1 
        kps4.forces[i] .= zeros(3)
    end
    KiteModels.init_springs(kps4)
    return pos, vel, posd, veld
end

function plot2d(x, z)
    plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    plot!(x, z, seriestype = :scatter)
end

KiteModels.set_depower_steering(kps4, 0.25, 0.0)

x0 = Float64[] 
z0 = Float64[]
if typeof(kps4) <: KPS4
    pos, vel, posd, veld = init2()
    kps4.alpha_depower = deg2rad(2.2095658807330962) # from one point simulation
    kps4.v_wind_gnd .= [9.0, 0.0, 0.0]
    height = 134.14733504839947                      # from one point simulation
    kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)    
    kps4.stiffness_factor = 0.1
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
     push!(x, kps4.pos[i][1])
     push!(z, kps4.pos[i][3])
end

println("kite distance: $(norm(kps4.pos[7]))")
println(KiteModels.spring_forces(kps4))
println("alpha_depower [deg]: $(rad2deg(kps4.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")

plot2d(x, z)

