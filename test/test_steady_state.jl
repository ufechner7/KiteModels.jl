using StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps
    const kps = KPS4(kcu)
end

const dt = 0.05

clear!(kps)
KiteModels.set_depower_steering!(kps, kps.set.depower_offset/100.0, 0.0)
kps.stiffness_factor = 0.5

@time KiteModels.find_steady_state!(kps, prn=true)

println("\nlift, drag    [N]  : $(KiteModels.lift_drag(kps))")
# println("\nSpring forces:")
# spring_forces(kps)
