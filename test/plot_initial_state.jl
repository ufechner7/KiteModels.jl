const FRONT_VIEW = false

clear(kps4)
KiteModels.set_depower_steering(kps4, kps4.set.depower_offset/100.0, 0.0)
height = sin(deg2rad(kps4.set.elevation)) * kps4.set.l_tether
kps4.v_wind .= kps4.v_wind_gnd * calc_wind_factor(kps4, height)
if typeof(kps4) <: KPS4
    kps4.stiffness_factor = 0.04
end

y0, yd0 = KiteModels.init(kps4)

find_steady_state(kps4, true)

if typeof(kps4) <: KPS4
    println("kite distance: $(norm(kps4.pos[end-2]))")
else
    println("kite distance: $(norm(kps4.pos[end]))")
end
println(KiteModels.spring_forces(kps4))
println("alpha_depower [deg]: $(rad2deg(kps4.alpha_depower))")
println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")

include("plot2d.jl")
plot2d(kps4.pos; zoom=false)
