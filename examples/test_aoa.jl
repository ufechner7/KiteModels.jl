# test case for the angle of attack calculation

using Printf
using KiteModels, KitePodModels, KiteUtils, LinearAlgebra, Rotations

set = deepcopy(load_settings("system.yaml"))
set.elevation = 70.0
set.alpha_zero = 0.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using ControlPlots
plt.close("all")

kcu::KCU = KCU(set)
kps4::KPS4 = KPS4(kcu)
KiteModels.init_springs!(kps4)
KiteModels.init_masses!(kps4)
pos, vel, acc = KiteModels.init_pos_vel_acc(kps4)

function calc_aoa(s::KPS4, pos, vel; alpha_depower=0.0, rel_steering=0.0, old=false)
    # pos_B, pos_C, pos_D: position of the kite particles B, C, and D
    # v_B,   v_C,   v_D:   velocity of the kite particles B, C, and D
    pos_B, pos_C, pos_D = pos[s.set.segments+3], pos[s.set.segments+4], pos[s.set.segments+5]
    v_B,   v_C,   v_D   = vel[s.set.segments+3], vel[s.set.segments+4], vel[s.set.segments+5]
    va_2,  va_3,  va_4  = s.v_wind - v_B, s.v_wind - v_C, s.v_wind - v_D
 
    pos_centre = 0.5 * (pos_C + pos_D)
    delta = pos_B - pos_centre
    z = -normalize(delta)
    if old
        y = normalize(pos_C - pos_D)
    else
        y = normalize(pos_D - pos_C)
    end

    x = y × z
    s.x .= x; s.y .= y; s.z .= z # save the kite reference frame in the state

    va_xz2 = va_2 - (va_2 ⋅ y) * y
    va_xy3 = va_3 - (va_3 ⋅ z) * z
    va_xy4 = va_4 - (va_4 ⋅ z) * z
    println("old: $old, x: $x, va_xz2: $va_xz2", " va_xy3: $va_xy3", " va_xy4: $va_xy4")
    if old
        alpha_2 = rad2deg(π - acos2(normalize(va_xz2) ⋅ x) - alpha_depower)     + s.set.alpha_zero
        alpha_3 = rad2deg(π - acos2(normalize(va_xy3) ⋅ x) - rel_steering * s.ks) + s.set.alpha_ztip
        alpha_4 = rad2deg(π - acos2(normalize(va_xy4) ⋅ x) + rel_steering * s.ks) + s.set.alpha_ztip
    else
        alpha_2 = rad2deg(π - acos2(normalize(va_xz2) ⋅ -x) - alpha_depower)     + s.set.alpha_zero
        alpha_3 = rad2deg(π - acos2(normalize(va_xy3) ⋅ -x) - rel_steering * s.ks) + s.set.alpha_ztip
        alpha_4 = rad2deg(π - acos2(normalize(va_xy4) ⋅ -x) + rel_steering * s.ks) + s.set.alpha_ztip
    end
    alpha_2, alpha_3, alpha_4
end


reltime=0.0
zoom=false
# plot2d(kps4.pos, reltime; zoom, xlim=(0,60), front=false, segments=set.segments)
vel[kps4.set.segments+3] .= [1.0, 2.0, 3.0]
vel[kps4.set.segments+4] .= [1.0, 2.0, 3.0]
vel[kps4.set.segments+5] .= [1.0, 2.0, 3.0]
pos[kps4.set.segments+3] .+= [0.1, 0.2, 0.3]
pos[kps4.set.segments+4] .+= [0.1, 0.2, 0.3]
pos[kps4.set.segments+5] .+= [0.1, 0.2, 0.3]
res=calc_aoa(kps4, pos, vel; old=true)
println("alpha_2: ", res[1], " alpha_3: ", res[2], " alpha_4: ", res[3])
res=calc_aoa(kps4, pos, vel; old=false)
println("alpha_2: ", res[1], " alpha_3: ", res[2], " alpha_4: ", res[3])