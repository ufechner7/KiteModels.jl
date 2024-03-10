using StaticArrays, LinearAlgebra, KiteUtils

using KiteModels, KitePodModels

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU(se())
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

res1 = zeros(SVector{SEGMENTS+1, KiteModels.KVec3})
res2 = deepcopy(res1)
if ! @isdefined res3
    const res3 = vcat(reduce(vcat, vcat(res1, res2)), zeros(2))
end

STEPS = 500
kps4.set.depower = 23.6
integrator = KiteModels.init_sim!(kps4; stiffness_factor=0.035, prn=false)
