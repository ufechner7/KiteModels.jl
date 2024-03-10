using KiteUtils, KiteModels, KitePodModels

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU(se())
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

kps4.set.depower = 23.6
integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
nothing
